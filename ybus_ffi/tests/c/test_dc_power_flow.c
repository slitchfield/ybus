/* The DC power flow: a solution's columns by the case's field codes, the solver, the errors,
 * and case14 against MATPOWER's committed results. */

#include <string.h>

#include "cases.h"
#include "check.h"

/* MATPOWER's results agree to within rounding: ybus_powerflow's own reference tests allow 1e-9 */
#define TOLERANCE 1e-9

/* case14's DC power flow against MATPOWER's rundcpf */
static void case14(void) {
    const char *path = matpower_case("case14");
    if (path == NULL) {
        return;
    }
    ybus_case *case14 = NULL;
    ybus_dc_solution *solution = NULL;
    CHECK(ybus_case_from_path(path, &case14) == YBUS_OK);
    CHECK(ybus_dc_power_flow(case14, NULL, &solution) == YBUS_OK);
    double va[14], pg[5], pf[20], expected[20];
    CHECK(ybus_dc_solution_get_bus(solution, YBUS_BUS_VA_DEG, va, 14) == YBUS_OK);
    CHECK(ybus_dc_solution_get_gen(solution, YBUS_GEN_PG_MW, pg, 5) == YBUS_OK);
    CHECK(ybus_dc_solution_get_branch(solution, YBUS_BRANCH_PF_MW, pf, 20) == YBUS_OK);
    CHECK(read_reference("case14", "dcpf", "bus", "va_deg", expected, 20) == 14);
    CHECK(all_close(va, expected, 14, TOLERANCE));
    CHECK(read_reference("case14", "dcpf", "gen", "pg_mw", expected, 20) == 5);
    CHECK(all_close(pg, expected, 5, TOLERANCE));
    CHECK(read_reference("case14", "dcpf", "branch", "pf_mw", expected, 20) == 20);
    CHECK(all_close(pf, expected, 20, TOLERANCE));
    ybus_dc_solution_free(solution);
    ybus_case_free(case14);
}

int main(void) {
    ybus_case *three = three_bus();

    /* The defaults (NULL options): the reference angle stays the case's, and the generators
     * cover the loads and the shunt, since DC power flow is lossless */
    ybus_dc_solution *solution = NULL;
    CHECK(ybus_dc_power_flow(three, NULL, &solution) == YBUS_OK);
    double va[3], pg[2], pf[3];
    CHECK(ybus_dc_solution_get_bus(solution, YBUS_BUS_VA_DEG, va, 3) == YBUS_OK);
    CHECK(ybus_dc_solution_get_gen(solution, YBUS_GEN_PG_MW, pg, 2) == YBUS_OK);
    CHECK(ybus_dc_solution_get_branch(solution, YBUS_BRANCH_PF_MW, pf, 3) == YBUS_OK);
    CHECK(va[0] == 0.0);
    double total = pg[0] + pg[1] - (0 + 50 + 80) - 1.5;
    CHECK(total < 1e-9 && total > -1e-9);
    ybus_dc_solver solver = 0;
    CHECK(ybus_dc_solution_solver(solution, &solver) == YBUS_OK);
    CHECK(solver == YBUS_DC_SOLVER_CHOLESKY);

    /* Only the solution's columns, and room for every row */
    CHECK(ybus_dc_solution_get_bus(solution, YBUS_BUS_VM_PU, va, 3) == YBUS_ERROR_INVALID_FIELD);
    CHECK(ybus_dc_solution_get_branch(solution, YBUS_BRANCH_PF_MW, pf, 2) ==
          YBUS_ERROR_BUFFER_TOO_SMALL);

    /* COLAMD factors by LU, to the same solution */
    ybus_dc_options *options = NULL;
    ybus_dc_solution *by_lu = NULL;
    CHECK(ybus_dc_options_new(&options) == YBUS_OK);
    CHECK(ybus_dc_options_set_ordering(options, YBUS_ORDERING_COLAMD) == YBUS_OK);
    CHECK(ybus_dc_power_flow(three, options, &by_lu) == YBUS_OK);
    CHECK(ybus_dc_solution_solver(by_lu, &solver) == YBUS_OK && solver == YBUS_DC_SOLVER_LU);
    double va_lu[3];
    CHECK(ybus_dc_solution_get_bus(by_lu, YBUS_BUS_VA_DEG, va_lu, 3) == YBUS_OK);
    CHECK(all_close(va_lu, va, 3, 1e-12));
    ybus_dc_solution_free(by_lu);
    ybus_dc_options_free(options);
    ybus_dc_solution_free(solution);
    ybus_dc_solution_free(NULL);

    /* With every bus PQ, nothing can be the reference */
    const int64_t pq[3] = {1, 1, 1};
    CHECK(ybus_case_set_bus_i64(three, YBUS_BUS_TYPE, pq, 3) == YBUS_OK);
    solution = (ybus_dc_solution *)three; /* any non-null value: it's set to NULL first */
    CHECK(ybus_dc_power_flow(three, NULL, &solution) == YBUS_ERROR_NO_REFERENCE_BUS);
    CHECK(solution == NULL);
    CHECK(ybus_dc_power_flow(NULL, NULL, &solution) == YBUS_ERROR_NULL_ARGUMENT);
    CHECK(strcmp(last_message(), "`case` is null") == 0);
    ybus_case_free(three);

    case14();
    return check_result();
}
