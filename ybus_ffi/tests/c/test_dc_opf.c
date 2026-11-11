/* The DC optimal power flow: a solution's columns by the case's field codes and its scalars,
 * a solve that stops short, and case14 against MATPOWER's committed rundcopf results. */

#include <math.h>
#include <string.h>

#include "cases.h"
#include "check.h"

/* MATPOWER's results agree to within rounding: powerflow's own reference tests allow 1e-9,
 * and 1e-12 relative for the objective */
#define TOLERANCE 1e-9
#define OBJECTIVE 1e-12

static int column_matches(ybus_status (*get)(const ybus_dc_opf_solution *, ybus_field,
                                             double *, size_t),
                          const ybus_dc_opf_solution *solution, ybus_field field,
                          const char *table, const char *column, size_t rows) {
    double values[64], expected[64];
    if (get(solution, field, values, rows) != YBUS_OK ||
        read_reference("case14", "dcopf", table, column, expected, 64) != rows) {
        return 0;
    }
    return all_close(values, expected, rows, TOLERANCE);
}

/* case14 against MATPOWER's rundcopf */
static void case14(void) {
    const char *path = matpower_case("case14");
    if (path == NULL) {
        return;
    }
    ybus_case *case14 = NULL;
    ybus_dc_opf_solution *solution = NULL;
    CHECK(ybus_case_from_path(path, &case14) == YBUS_OK);
    CHECK(ybus_dc_opf(case14, NULL, &solution) == YBUS_OK);
    CHECK(column_matches(ybus_dc_opf_solution_get_bus, solution, YBUS_BUS_VA_DEG, "bus",
                         "va_deg", 14));
    CHECK(column_matches(ybus_dc_opf_solution_get_bus, solution, YBUS_BUS_LAM_P, "bus", "lam_p",
                         14));
    CHECK(column_matches(ybus_dc_opf_solution_get_gen, solution, YBUS_GEN_PG_MW, "gen", "pg_mw",
                         5));
    CHECK(column_matches(ybus_dc_opf_solution_get_gen, solution, YBUS_GEN_MU_PMIN, "gen",
                         "mu_pmin", 5));
    CHECK(column_matches(ybus_dc_opf_solution_get_gen, solution, YBUS_GEN_MU_PMAX, "gen",
                         "mu_pmax", 5));
    CHECK(column_matches(ybus_dc_opf_solution_get_branch, solution, YBUS_BRANCH_PF_MW, "branch",
                         "pf_mw", 20));
    CHECK(column_matches(ybus_dc_opf_solution_get_branch, solution, YBUS_BRANCH_MU_SF, "branch",
                         "mu_sf", 20));
    CHECK(column_matches(ybus_dc_opf_solution_get_branch, solution, YBUS_BRANCH_MU_ST, "branch",
                         "mu_st", 20));
    CHECK(column_matches(ybus_dc_opf_solution_get_branch, solution, YBUS_BRANCH_MU_ANGMIN,
                         "branch", "mu_angmin", 20));
    CHECK(column_matches(ybus_dc_opf_solution_get_branch, solution, YBUS_BRANCH_MU_ANGMAX,
                         "branch", "mu_angmax", 20));
    double objective = 0, expected[1];
    size_t iterations = 0;
    ybus_opf_status status = 9;
    CHECK(ybus_dc_opf_solution_objective(solution, &objective) == YBUS_OK);
    CHECK(read_reference("case14", "dcopf", "summary", "objective", expected, 1) == 1);
    CHECK(fabs(objective - expected[0]) <= OBJECTIVE * fabs(expected[0]));
    CHECK(ybus_dc_opf_solution_iterations(solution, &iterations) == YBUS_OK);
    CHECK(read_reference("case14", "dcopf", "summary", "iterations", expected, 1) == 1);
    CHECK((double)iterations == expected[0]);
    CHECK(ybus_dc_opf_solution_status(solution, &status) == YBUS_OK);
    CHECK(read_reference("case14", "dcopf", "summary", "success", expected, 1) == 1);
    CHECK(status == YBUS_OPF_CONVERGED && expected[0] == 1);
    ybus_dc_opf_solution_free(solution);
    ybus_case_free(case14);
}

int main(void) {
    ybus_case *three = three_bus();

    /* MATPOWER's defaults (NULL options) converge; the generation covers the load and the
     * shunt, and each bus's price is set */
    ybus_dc_opf_solution *solution = NULL;
    CHECK(ybus_dc_opf(three, NULL, &solution) == YBUS_OK);
    ybus_opf_status status = 9;
    CHECK(ybus_dc_opf_solution_status(solution, &status) == YBUS_OK);
    CHECK(status == YBUS_OPF_CONVERGED);
    double pg[2], lam[3], objective = 0;
    CHECK(ybus_dc_opf_solution_get_gen(solution, YBUS_GEN_PG_MW, pg, 2) == YBUS_OK);
    CHECK(fabs(pg[0] + pg[1] - 131.5) < 1e-6);
    CHECK(ybus_dc_opf_solution_get_bus(solution, YBUS_BUS_LAM_P, lam, 3) == YBUS_OK);
    CHECK(lam[0] > 0);
    CHECK(ybus_dc_opf_solution_objective(solution, &objective) == YBUS_OK && objective > 0);
    CHECK(ybus_dc_opf_solution_get_bus(solution, YBUS_BUS_VM_PU, lam, 3) ==
          YBUS_ERROR_INVALID_FIELD);
    ybus_dc_opf_solution_free(solution);
    ybus_dc_opf_solution_free(NULL);

    /* One iteration isn't enough: a solution still, whose status says so */
    ybus_dc_opf_options *options = NULL;
    size_t iterations = 0;
    CHECK(ybus_dc_opf_options_new(YBUS_OPF_MIPS, &options) == YBUS_OK);
    CHECK(ybus_dc_opf_options_set_max_iterations(options, 1) == YBUS_OK);
    CHECK(ybus_dc_opf(three, options, &solution) == YBUS_OK);
    CHECK(ybus_dc_opf_solution_status(solution, &status) == YBUS_OK);
    CHECK(status == YBUS_OPF_MAX_ITERATIONS);
    CHECK(ybus_dc_opf_solution_iterations(solution, &iterations) == YBUS_OK && iterations == 1);
    ybus_dc_opf_solution_free(solution);
    ybus_dc_opf_options_free(options);

    CHECK(ybus_dc_opf(NULL, NULL, &solution) == YBUS_ERROR_NULL_ARGUMENT);
    CHECK(strcmp(last_message(), "`case` is null") == 0 && solution == NULL);
    ybus_case_free(three);

    case14();
    return check_result();
}
