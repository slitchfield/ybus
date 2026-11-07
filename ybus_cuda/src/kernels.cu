// ybus_cuda's kernels, compiled at run time by NVRTC. Every block is column-major: column k of
// an n-row block starts at k * n. A row index of -1 means "none".

// Sets each column of `block` (zeroed beforehand) to e_plus - e_minus.
extern "C" __global__ void scatter_pairs(double *block, int n, int width, const int *plus,
                                         const int *minus) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= width) {
        return;
    }
    if (plus[k] >= 0) {
        block[(long)k * n + plus[k]] = 1.0;
    }
    if (minus[k] >= 0) {
        block[(long)k * n + minus[k]] = -1.0;
    }
}

// Each column's rank-one update of `base`, given y = A^-1 c in that column of `y`, with
// c = e_plus - e_minus (read at the solution's rows `plus` and `minus`):
//   x = base + alpha y,  w = inverse_weight - c'y,  z = c'x / w,  out = x + z y,
// as the CPU's Sherman-Morrison-Woodbury update orders it for one column. A column with w = 0 is
// flagged in `singular` and gets x.
extern "C" __global__ void pair_update(const double *y, double *out, int n, int width,
                                       const double *base, const int *plus, const int *minus,
                                       const double *alpha, const double *inverse_weight,
                                       int *singular) {
    long index = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (long)n * width) {
        return;
    }
    int k = (int)(index / n);
    int i = (int)(index % n);
    const double *column = y + (long)k * n;
    int p = plus[k];
    int m = minus[k];
    double a = alpha[k];
    double cy = 0.0;
    double cx = 0.0;
    if (p >= 0) {
        cy += column[p];
        cx += base[p] + a * column[p];
    }
    if (m >= 0) {
        cy -= column[m];
        cx -= base[m] + a * column[m];
    }
    double w = -cy + inverse_weight[k];
    double x = base[i] + a * column[i];
    if (w == 0.0) {
        if (i == 0) {
            singular[k] = 1;
        }
        out[index] = x;
        return;
    }
    out[index] = x + (cx / w) * column[i];
}

// Sets each column of `block`'s entry at `rows[k]` (if not -1) to zero.
extern "C" __global__ void zero_entries(double *block, int n, int width, const int *rows) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < width && rows[k] >= 0) {
        block[(long)k * n + rows[k]] = 0.0;
    }
}

// Sets every column of `block` to `offset`.
extern "C" __global__ void broadcast(double *block, int n, int width, const double *offset) {
    long index = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (index < (long)n * width) {
        block[index] = offset[index % n];
    }
}

// One thread block per column: the largest |value| * weight[row] and its row, the smallest row on a tie, or -1 if every product is zero; and whether any value
// isn't finite. 256 threads per block.
extern "C" __global__ void column_max(const double *block, int n, int width,
                                      const double *weight, double *best, int *best_row,
                                      int *not_finite) {
    __shared__ double value[256];
    __shared__ int row[256];
    __shared__ int bad[256];
    int k = blockIdx.x;
    if (k >= width) {
        return;
    }
    const double *column = block + (long)k * n;
    double my_value = 0.0;
    int my_row = -1;
    int my_bad = 0;
    for (int i = threadIdx.x; i < n; i += blockDim.x) {
        double v = column[i];
        if (!isfinite(v)) {
            my_bad = 1;
            continue;
        }
        double product = fabs(v) * weight[i];
        // Rows are visited in increasing order, so a strictly larger product wins
        if (product > my_value) {
            my_value = product;
            my_row = i;
        }
    }
    value[threadIdx.x] = my_value;
    row[threadIdx.x] = my_row;
    bad[threadIdx.x] = my_bad;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
        if (threadIdx.x < stride) {
            int other = threadIdx.x + stride;
            double v = value[other];
            int r = row[other];
            if (r >= 0 && (v > value[threadIdx.x] ||
                           (v == value[threadIdx.x] && (row[threadIdx.x] < 0 || r < row[threadIdx.x])))) {
                value[threadIdx.x] = v;
                row[threadIdx.x] = r;
            }
            bad[threadIdx.x] |= bad[other];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        best[k] = value[0];
        best_row[k] = row[0];
        not_finite[k] = bad[0];
    }
}
