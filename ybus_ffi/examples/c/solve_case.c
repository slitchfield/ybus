/* solve_case: loads a MATPOWER case and runs each of ybus's solves on it, the DC and AC power
 * flows and the DC and AC optimal power flows, printing a few results of each.
 *
 * A first look at the C interface: each handle is made by a function and freed by its own
 * _free, every call returns a ybus_status, and a solution's columns are read whole into
 * arrays, by the same field codes as the case's columns (MATPOWER's column numbers).
 *
 * Build and run, from the repository root:
 *   cargo build -p ybus_ffi --release
 *   cc -I ybus_ffi/include ybus_ffi/examples/c/solve_case.c target/release/libybus.a \
 *       -lgcc_s -lutil -lrt -lpthread -lm -ldl -lc -o solve_case
 *   ./solve_case testdata/matpower/matpower8.1/data/case14.m
 */

#include <stdio.h>
#include <stdlib.h>

#include "ybus.h"

/* Prints why `what` failed: its status, and the message the library kept for this thread */
static void report(const char *what, ybus_status status) {
    char message[512];
    size_t needed = ybus_last_error_message(message, sizeof message);
    fprintf(stderr, "%s failed (status %d): %s\n", what, (int)status,
            needed > 0 && needed <= sizeof message ? message : "(no message)");
}

/* The smallest and largest of n values */
static void range(const double *values, size_t n, double *low, double *high) {
    *low = *high = values[0];
    for (size_t i = 1; i < n; i++) {
        if (values[i] < *low) {
            *low = values[i];
        }
        if (values[i] > *high) {
            *high = values[i];
        }
    }
}

/* How an OPF's solver ended, in words */
static const char *status_name(ybus_opf_status status) {
    switch (status) {
    case YBUS_OPF_CONVERGED:
        return "converged";
    case YBUS_OPF_MAX_ITERATIONS:
        return "stopped at its iteration limit";
    case YBUS_OPF_NUMERICAL_FAILURE:
        return "failed numerically";
    default:
        return "ended in an unknown way";
    }
}

/* The DC power flow, with default options: its range of bus angles */
static int dc_power_flow(const ybus_case *grid, double *buses, size_t bus_count) {
    ybus_dc_solution *solution = NULL;
    ybus_status status = ybus_dc_power_flow(grid, NULL, &solution);
    if (status != YBUS_OK) {
        report("The DC power flow", status);
        return 1;
    }
    status = ybus_dc_solution_get_bus(solution, YBUS_BUS_VA_DEG, buses, bus_count);
    if (status == YBUS_OK) {
        double low, high;
        range(buses, bus_count, &low, &high);
        printf("DC power flow: bus angles from %.3f to %.3f degrees\n", low, high);
    } else {
        report("Reading the DC angles", status);
    }
    ybus_dc_solution_free(solution);
    return status != YBUS_OK;
}

/* The AC power flow, by Newton's method: its iterations and range of voltage magnitudes */
static int ac_power_flow(const ybus_case *grid, double *buses, size_t bus_count) {
    ybus_ac_solution *solution = NULL;
    ybus_status status = ybus_ac_power_flow(grid, NULL, &solution);
    if (status != YBUS_OK) {
        report("The AC power flow", status);
        return 1;
    }
    size_t iterations = 0;
    status = ybus_ac_solution_iterations(solution, &iterations);
    if (status == YBUS_OK) {
        status = ybus_ac_solution_get_bus(solution, YBUS_BUS_VM_PU, buses, bus_count);
    }
    if (status == YBUS_OK) {
        double low, high;
        range(buses, bus_count, &low, &high);
        printf("AC power flow: %zu Newton iterations, voltages from %.4f to %.4f p.u.\n",
               iterations, low, high);
    } else {
        report("Reading the AC solution", status);
    }
    ybus_ac_solution_free(solution);
    return status != YBUS_OK;
}

/* The DC optimal power flow, with MATPOWER's defaults: how it ended, its cost and its highest
 * real power price */
static int dc_opf(const ybus_case *grid, double *buses, size_t bus_count) {
    ybus_dc_opf_solution *solution = NULL;
    ybus_status status = ybus_dc_opf(grid, NULL, &solution);
    if (status != YBUS_OK) {
        report("The DC OPF", status);
        return 1;
    }
    ybus_opf_status ended = 0;
    double objective = 0, low, high;
    size_t iterations = 0;
    if ((status = ybus_dc_opf_solution_status(solution, &ended)) == YBUS_OK &&
        (status = ybus_dc_opf_solution_objective(solution, &objective)) == YBUS_OK &&
        (status = ybus_dc_opf_solution_iterations(solution, &iterations)) == YBUS_OK &&
        (status = ybus_dc_opf_solution_get_bus(solution, YBUS_BUS_LAM_P, buses, bus_count)) ==
            YBUS_OK) {
        range(buses, bus_count, &low, &high);
        printf("DC OPF: %s in %zu iterations, cost %.2f $/hr, prices up to %.2f $/MWh\n",
               status_name(ended), iterations, objective, high);
    } else {
        report("Reading the DC OPF solution", status);
    }
    ybus_dc_opf_solution_free(solution);
    return status != YBUS_OK || ended != YBUS_OPF_CONVERGED;
}

/* The AC optimal power flow, with MATPOWER's defaults: the same, with reactive power priced
 * too */
static int ac_opf(const ybus_case *grid, double *buses, size_t bus_count) {
    ybus_ac_opf_solution *solution = NULL;
    ybus_status status = ybus_ac_opf(grid, NULL, &solution);
    if (status != YBUS_OK) {
        report("The AC OPF", status);
        return 1;
    }
    ybus_opf_status ended = 0;
    double objective = 0, low, high;
    size_t iterations = 0;
    if ((status = ybus_ac_opf_solution_status(solution, &ended)) == YBUS_OK &&
        (status = ybus_ac_opf_solution_objective(solution, &objective)) == YBUS_OK &&
        (status = ybus_ac_opf_solution_iterations(solution, &iterations)) == YBUS_OK &&
        (status = ybus_ac_opf_solution_get_bus(solution, YBUS_BUS_LAM_P, buses, bus_count)) ==
            YBUS_OK) {
        range(buses, bus_count, &low, &high);
        printf("AC OPF: %s in %zu iterations, cost %.2f $/hr, prices up to %.2f $/MWh\n",
               status_name(ended), iterations, objective, high);
    } else {
        report("Reading the AC OPF solution", status);
    }
    ybus_ac_opf_solution_free(solution);
    return status != YBUS_OK || ended != YBUS_OPF_CONVERGED;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s <case.m>\n", argv[0]);
        return 2;
    }
    ybus_case *grid = NULL;
    ybus_status status = ybus_case_from_path(argv[1], &grid);
    if (status != YBUS_OK) {
        report("Loading the case", status);
        return 1;
    }
    size_t buses = 0, generators = 0, branches = 0;
    ybus_case_bus_count(grid, &buses);
    ybus_case_gen_count(grid, &generators);
    ybus_case_branch_count(grid, &branches);
    printf("%s: %zu buses, %zu generators, %zu branches\n", argv[1], buses, generators,
           branches);

    /* One buffer, a value per bus, for every bus column read below */
    double *column = malloc(buses * sizeof *column);
    if (column == NULL) {
        ybus_case_free(grid);
        return 1;
    }
    /* Each solve reports its own failure and the others still run */
    int failures = dc_power_flow(grid, column, buses) + ac_power_flow(grid, column, buses) +
                   dc_opf(grid, column, buses) + ac_opf(grid, column, buses);
    free(column);
    ybus_case_free(grid);
    return failures > 0;
}
