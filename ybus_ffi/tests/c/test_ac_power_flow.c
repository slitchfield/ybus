/* The AC power flow: a solution's columns by the case's field codes and its scalars, the
 * errors, and case14 against MATPOWER's committed results, without and with reactive limits. */

#include <string.h>

#include "cases.h"
#include "check.h"

/* MATPOWER's results agree to within rounding: powerflow's own reference tests allow 1e-9 */
#define TOLERANCE 1e-9

/* One column of a solution against the reference's column of that name */
static int column_matches(ybus_status (*get)(const ybus_ac_solution *, ybus_field, double *,
                                             size_t),
                          const ybus_ac_solution *solution, ybus_field field, const char *run,
                          const char *table, const char *column, size_t rows) {
    double values[64], expected[64];
    if (get(solution, field, values, rows) != YBUS_OK ||
        read_reference("case14", run, table, column, expected, 64) != rows) {
        return 0;
    }
    return all_close(values, expected, rows, TOLERANCE);
}

/* case14 against MATPOWER's runpf: Newton, then with reactive limits enforced all at once */
static void case14(void) {
    const char *path = matpower_case("case14");
    if (path == NULL) {
        return;
    }
    ybus_case *case14 = NULL;
    ybus_ac_solution *solution = NULL;
    CHECK(ybus_case_from_path(path, &case14) == YBUS_OK);
    CHECK(ybus_ac_power_flow(case14, NULL, &solution) == YBUS_OK);
    CHECK(column_matches(ybus_ac_solution_get_bus, solution, YBUS_BUS_VM_PU, "acpf", "bus",
                         "vm_pu", 14));
    CHECK(column_matches(ybus_ac_solution_get_bus, solution, YBUS_BUS_VA_DEG, "acpf", "bus",
                         "va_deg", 14));
    CHECK(column_matches(ybus_ac_solution_get_gen, solution, YBUS_GEN_PG_MW, "acpf", "gen",
                         "pg_mw", 5));
    CHECK(column_matches(ybus_ac_solution_get_gen, solution, YBUS_GEN_QG_MVAR, "acpf", "gen",
                         "qg_mvar", 5));
    CHECK(column_matches(ybus_ac_solution_get_branch, solution, YBUS_BRANCH_PF_MW, "acpf",
                         "branch", "pf_mw", 20));
    CHECK(column_matches(ybus_ac_solution_get_branch, solution, YBUS_BRANCH_QF_MVAR, "acpf",
                         "branch", "qf_mvar", 20));
    CHECK(column_matches(ybus_ac_solution_get_branch, solution, YBUS_BRANCH_PT_MW, "acpf",
                         "branch", "pt_mw", 20));
    CHECK(column_matches(ybus_ac_solution_get_branch, solution, YBUS_BRANCH_QT_MVAR, "acpf",
                         "branch", "qt_mvar", 20));
    double expected[1];
    size_t iterations = 0;
    CHECK(ybus_ac_solution_iterations(solution, &iterations) == YBUS_OK);
    CHECK(read_reference("case14", "acpf", "summary", "iterations", expected, 1) == 1);
    CHECK((double)iterations == expected[0]);
    ybus_ac_solution_free(solution);

    /* Reactive limits, all at once: MATPOWER numbers the limited generators from 1 */
    ybus_ac_options *options = NULL;
    CHECK(ybus_ac_options_new(YBUS_AC_NEWTON, &options) == YBUS_OK);
    CHECK(ybus_ac_options_set_reactive_limits(options, YBUS_REACTIVE_LIMITS_ALL_AT_ONCE) ==
          YBUS_OK);
    CHECK(ybus_ac_power_flow(case14, options, &solution) == YBUS_OK);
    CHECK(column_matches(ybus_ac_solution_get_gen, solution, YBUS_GEN_QG_MVAR, "qlim1", "gen",
                         "qg_mvar", 5));
    size_t count = 0, limited[5];
    CHECK(ybus_ac_solution_limited_generators(solution, NULL, 0, &count) == YBUS_OK);
    CHECK(read_reference("case14", "qlim1", "limited", "row", expected, 1) == count);
    CHECK(ybus_ac_solution_limited_generators(solution, limited, 5, &count) == YBUS_OK);
    CHECK(count == 1 && (double)(limited[0] + 1) == expected[0]);
    CHECK(ybus_ac_solution_iterations(solution, &iterations) == YBUS_OK);
    CHECK(read_reference("case14", "qlim1", "summary", "iterations", expected, 1) == 1);
    CHECK((double)iterations == expected[0]);
    ybus_ac_solution_free(solution);
    ybus_ac_options_free(options);
    ybus_case_free(case14);
}

int main(void) {
    ybus_case *three = three_bus();

    /* Newton's defaults (NULL options) converge, within the tolerance */
    ybus_ac_solution *solution = NULL;
    CHECK(ybus_ac_power_flow(three, NULL, &solution) == YBUS_OK);
    double vm[3], mismatch = 1;
    size_t iterations = 0, factorizations = 0, count = 9;
    CHECK(ybus_ac_solution_get_bus(solution, YBUS_BUS_VM_PU, vm, 3) == YBUS_OK);
    CHECK(vm[0] == 1.02 && vm[1] == 1.0);
    CHECK(ybus_ac_solution_iterations(solution, &iterations) == YBUS_OK && iterations > 0);
    CHECK(ybus_ac_solution_mismatch_pu(solution, &mismatch) == YBUS_OK && mismatch < 1e-8);
    CHECK(ybus_ac_solution_full_factorizations(solution, &factorizations) == YBUS_OK);
    CHECK(factorizations == 1);
    CHECK(ybus_ac_solution_limited_generators(solution, NULL, 0, &count) == YBUS_OK);
    CHECK(count == 0);
    CHECK(ybus_ac_solution_get_bus(solution, YBUS_BUS_PD_MW, vm, 3) == YBUS_ERROR_INVALID_FIELD);
    ybus_ac_solution_free(solution);
    ybus_ac_solution_free(NULL);

    /* No iterations allowed: it can't converge */
    ybus_ac_options *options = NULL;
    CHECK(ybus_ac_options_new(YBUS_AC_NEWTON, &options) == YBUS_OK);
    CHECK(ybus_ac_options_set_max_iterations(options, 0) == YBUS_OK);
    solution = (ybus_ac_solution *)three; /* any non-null value: it's set to NULL first */
    CHECK(ybus_ac_power_flow(three, options, &solution) == YBUS_ERROR_NO_CONVERGENCE);
    CHECK(solution == NULL && strlen(last_message()) > 0);
    ybus_ac_options_free(options);
    ybus_case_free(three);

    case14();
    return check_result();
}
