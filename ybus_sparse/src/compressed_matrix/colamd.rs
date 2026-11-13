//! COLAMD: a column ordering for LU of a matrix whose pattern needn't be symmetric (Davis,
//! Gilbert, Larimore and Ng, 2004).
//!
//! The columns of `A Q` are ordered so the Cholesky factor of `(A Q)ᵀ (A Q)` stays sparse, which
//! bounds the fill of `A Q`'s LU factors under any row pivoting, without forming `AᵀA`. The
//! method works on `A`'s row–column structure:
//!
//! - **Rows as elements:** each row is the set of columns it joins. Eliminating a column `c`
//!   merges every row holding `c` into one new *pivot row*, the union of their columns (less
//!   `c`): the pattern of `AᵀA`'s elimination stays implicit in these unions.
//! - **Approximate scores:** a column's score approximates its external degree in `AᵀA`, the
//!   size of the pivot row it would create: for column `j`, the pivot row's other columns plus,
//!   for each other row of `j`, the columns it holds outside the pivot row (`|r \ p|`), counted
//!   with one pass over the pivot row's columns. The column of smallest score is eliminated
//!   next.
//! - **Aggressive absorption:** a row with no columns outside the pivot row is absorbed into it.
//! - **Supercolumns:** columns of the pivot row whose rows are then the same are merged and
//!   ordered together, so each is scored once.
//! - **Dense rows and columns** (more than `max(16, 10√n)` entries, as this crate's AMD treats
//!   dense rows) are set aside: dense rows are ignored, and dense columns, with columns left
//!   with no rows, are ordered last.
//!
//! Written from the paper's description; no code from other implementations was used.

use std::collections::HashMap;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Status {
    Live,
    Eliminated,
    Merged,
    SetAside,
}

/// Degree lists: columns bucketed by score, each bucket a doubly linked list.
struct Buckets {
    head: Vec<Option<usize>>,
    next: Vec<Option<usize>>,
    prev: Vec<Option<usize>>,
    score: Vec<usize>,
    lowest: usize,
}

impl Buckets {
    fn new(n_cols: usize) -> Buckets {
        Buckets {
            head: vec![None; n_cols + 1],
            next: vec![None; n_cols],
            prev: vec![None; n_cols],
            score: vec![0; n_cols],
            lowest: n_cols,
        }
    }

    fn insert(&mut self, col: usize, score: usize) {
        self.score[col] = score;
        self.prev[col] = None;
        self.next[col] = self.head[score];
        if let Some(first) = self.head[score] {
            self.prev[first] = Some(col);
        }
        self.head[score] = Some(col);
        self.lowest = self.lowest.min(score);
    }

    fn remove(&mut self, col: usize) {
        let score = self.score[col];
        match self.prev[col] {
            Some(before) => self.next[before] = self.next[col],
            None => self.head[score] = self.next[col],
        }
        if let Some(after) = self.next[col] {
            self.prev[after] = self.prev[col];
        }
    }

    /// The column of lowest score (the most recently inserted, on a tie), if any.
    fn lowest(&mut self) -> Option<usize> {
        while self.lowest < self.head.len() {
            if let Some(col) = self.head[self.lowest] {
                return Some(col);
            }
            self.lowest += 1;
        }
        None
    }
}

/// The COLAMD column ordering of an `n_rows × n_cols` pattern given by its `(row, col)`
/// entries (duplicates allowed): `order[k]` is the column placed `k`-th.
pub(crate) fn column_approximate_minimum_degree(
    n_rows: usize,
    n_cols: usize,
    entries: impl IntoIterator<Item = (usize, usize)>,
) -> Vec<usize> {
    let mut col_rows: Vec<Vec<usize>> = vec![vec![]; n_cols];
    for (row, col) in entries {
        col_rows[col].push(row);
    }
    for rows in &mut col_rows {
        rows.sort_unstable();
        rows.dedup();
    }

    // Dense columns are set aside; then dense rows (counting the remaining columns) are
    //   ignored; columns left with no rows are ordered last, before the dense ones
    let threshold = |n: usize| 16.max((10.0 * (n as f64).sqrt()) as usize);
    let dense_col = threshold(n_rows.min(n_cols));
    let dense_row = threshold(n_cols);
    let mut status = vec![Status::Live; n_cols];
    let mut dense_cols = vec![];
    for (col, rows) in col_rows.iter().enumerate() {
        if rows.len() > dense_col {
            status[col] = Status::SetAside;
            dense_cols.push(col);
        }
    }
    let mut row_cols: Vec<Vec<usize>> = vec![vec![]; n_rows];
    for (col, rows) in col_rows.iter().enumerate() {
        if status[col] == Status::Live {
            for &row in rows {
                row_cols[row].push(col);
            }
        }
    }
    let mut row_alive: Vec<bool> = row_cols
        .iter()
        .map(|cols| !cols.is_empty() && cols.len() <= dense_row)
        .collect();
    let mut null_cols = vec![];
    for (col, rows) in col_rows.iter_mut().enumerate() {
        rows.retain(|&row| row_alive[row]);
        if status[col] == Status::Live && rows.is_empty() {
            status[col] = Status::SetAside;
            null_cols.push(col);
        }
    }

    let mut weight = vec![1; n_cols];
    let mut members: Vec<Vec<usize>> = vec![vec![]; n_cols];
    let mut row_size: Vec<usize> = row_cols.iter().map(Vec::len).collect();
    let mut live_weight = status.iter().filter(|&&s| s == Status::Live).count();
    let mut buckets = Buckets::new(n_cols);
    for col in 0..n_cols {
        if status[col] == Status::Live {
            let score: usize = col_rows[col].iter().map(|&row| row_size[row] - 1).sum();
            buckets.insert(col, score.min(live_weight - 1));
        }
    }

    let mut order = Vec::with_capacity(n_cols);
    let mut col_mark = vec![0usize; n_cols];
    let mut row_mark: Vec<usize> = vec![0; n_rows];
    let mut outside: Vec<isize> = vec![0; n_rows];
    let mut stamp = 0;
    while let Some(pivot) = buckets.lowest() {
        buckets.remove(pivot);
        status[pivot] = Status::Eliminated;
        live_weight -= weight[pivot];
        order.push(pivot);
        order.extend_from_slice(&members[pivot]);

        // The pivot row: the union of the pivot column's rows, which it absorbs
        stamp += 1;
        let mut pivot_cols = vec![];
        for &row in &col_rows[pivot] {
            if !row_alive[row] {
                continue;
            }
            for &col in &row_cols[row] {
                if status[col] == Status::Live && col_mark[col] != stamp {
                    col_mark[col] = stamp;
                    pivot_cols.push(col);
                }
            }
            row_alive[row] = false;
        }
        let p = row_cols.len();
        row_size.push(pivot_cols.iter().map(|&col| weight[col]).sum());
        row_cols.push(pivot_cols.clone());
        row_alive.push(true);
        row_mark.push(0);
        outside.push(0);
        for &col in &pivot_cols {
            buckets.remove(col);
            col_rows[col].retain(|&row| row_alive[row]);
            col_rows[col].push(p);
        }

        // |r \ p| for each row meeting the pivot row: its size less its columns in the pivot
        //   row, weighted
        for &col in &pivot_cols {
            for &row in &col_rows[col] {
                if row == p {
                    continue;
                }
                if row_mark[row] != stamp {
                    row_mark[row] = stamp;
                    outside[row] = row_size[row] as isize;
                }
                outside[row] -= weight[col] as isize;
            }
        }
        // Aggressive absorption: a row inside the pivot row is absorbed into it
        for &col in &pivot_cols {
            for &row in &col_rows[col] {
                if row != p && row_mark[row] == stamp && outside[row] == 0 {
                    row_alive[row] = false;
                }
            }
        }
        let mut scores = Vec::with_capacity(pivot_cols.len());
        for &col in &pivot_cols {
            col_rows[col].retain(|&row| row_alive[row]);
            let external: isize = col_rows[col]
                .iter()
                .filter(|&&row| row != p)
                .map(|&row| outside[row])
                .sum();
            scores.push(row_size[p] as isize - weight[col] as isize + external);
        }

        // Supercolumns: columns of the pivot row with the same rows, merged into the first
        let mut representatives: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (index, &col) in pivot_cols.iter().enumerate() {
            let mut rows = col_rows[col].clone();
            rows.sort_unstable();
            col_rows[col] = rows;
            let key = (
                col_rows[col].len(),
                col_rows[col].iter().sum::<usize>() % (n_cols + 1),
            );
            let candidates = representatives.entry(key).or_default();
            match candidates
                .iter()
                .find(|&&other| col_rows[pivot_cols[other]] == col_rows[col])
            {
                Some(&other) => {
                    let into = pivot_cols[other];
                    weight[into] += weight[col];
                    let mut moved = vec![col];
                    moved.append(&mut members[col]);
                    members[into].extend(moved);
                    status[col] = Status::Merged;
                    // The merged column is no longer outside the representative's rows
                    scores[other] -= weight[col] as isize;
                }
                None => candidates.push(index),
            }
        }
        for (&col, &score) in pivot_cols.iter().zip(&scores) {
            if status[col] == Status::Live {
                let cap = live_weight - weight[col];
                buckets.insert(col, (score.max(0) as usize).min(cap));
            }
        }
    }
    order.extend(null_cols);
    order.extend(dense_cols);
    order
}
