/* The AC optimal power flow: a generator it can't use refused, a solution's columns by the
 * case's field codes and its scalars, and case14 against MATPOWER's committed runopf results. */

#include <math.h>
#include <string.h>

#include "cases.h"
#include "check.h"

/* MATPOWER's results agree to within rounding: ybus_powerflow's own reference tests allow 1e-8 for
 * the AC OPF's columns, its multipliers included, and 1e-12 relative for the objective */
#define TOLERANCE 1e-8
#define OBJECTIVE 1e-12

static int column_matches(ybus_status (*get)(const ybus_ac_opf_solution *, ybus_field,
                                             double *, size_t),
                          const ybus_ac_opf_solution *solution, ybus_field field,
                          const char *table, const char *column, size_t rows) {
    double values[64], expected[64];
    if (get(solution, field, values, rows) != YBUS_OK ||
        read_reference("case14", "opf", table, column, expected, 64) != rows) {
        return 0;
    }
    return all_close(values, expected, rows, TOLERANCE);
}

/* case14 against MATPOWER's runopf, every column it reports */
static void case14(void) {
    const char *path = matpower_case("case14");
    if (path == NULL) {
        return;
    }
    ybus_case *case14 = NULL;
    ybus_ac_opf_solution *s = NULL;
    CHECK(ybus_case_from_path(path, &case14) == YBUS_OK);
    CHECK(ybus_ac_opf(case14, NULL, &s) == YBUS_OK);
    const struct {
        ybus_field field;
        const char *column;
    } buses[] = {{YBUS_BUS_VM_PU, "vm_pu"},     {YBUS_BUS_VA_DEG, "va_deg"},
                 {YBUS_BUS_LAM_P, "lam_p"},     {YBUS_BUS_LAM_Q, "lam_q"},
                 {YBUS_BUS_MU_VMIN, "mu_vmin"}, {YBUS_BUS_MU_VMAX, "mu_vmax"}},
      gens[] = {{YBUS_GEN_PG_MW, "pg_mw"},     {YBUS_GEN_QG_MVAR, "qg_mvar"},
                {YBUS_GEN_MU_PMIN, "mu_pmin"}, {YBUS_GEN_MU_PMAX, "mu_pmax"},
                {YBUS_GEN_MU_QMIN, "mu_qmin"}, {YBUS_GEN_MU_QMAX, "mu_qmax"}},
      branches[] = {{YBUS_BRANCH_PF_MW, "pf_mw"},         {YBUS_BRANCH_QF_MVAR, "qf_mvar"},
                    {YBUS_BRANCH_PT_MW, "pt_mw"},         {YBUS_BRANCH_QT_MVAR, "qt_mvar"},
                    {YBUS_BRANCH_MU_SF, "mu_sf"},         {YBUS_BRANCH_MU_ST, "mu_st"},
                    {YBUS_BRANCH_MU_ANGMIN, "mu_angmin"}, {YBUS_BRANCH_MU_ANGMAX, "mu_angmax"}};
    for (size_t i = 0; i < 6; i++) {
        CHECK(column_matches(ybus_ac_opf_solution_get_bus, s, buses[i].field, "bus",
                             buses[i].column, 14));
        CHECK(column_matches(ybus_ac_opf_solution_get_gen, s, gens[i].field, "gen",
                             gens[i].column, 5));
    }
    for (size_t i = 0; i < 8; i++) {
        CHECK(column_matches(ybus_ac_opf_solution_get_branch, s, branches[i].field, "branch",
                             branches[i].column, 20));
    }
    double objective = 0, expected[1];
    size_t iterations = 0;
    ybus_opf_status status = 9;
    CHECK(ybus_ac_opf_solution_objective(s, &objective) == YBUS_OK);
    CHECK(read_reference("case14", "opf", "summary", "objective", expected, 1) == 1);
    CHECK(fabs(objective - expected[0]) <= OBJECTIVE * fabs(expected[0]));
    CHECK(ybus_ac_opf_solution_iterations(s, &iterations) == YBUS_OK);
    CHECK(read_reference("case14", "opf", "summary", "iterations", expected, 1) == 1);
    CHECK((double)iterations == expected[0]);
    CHECK(ybus_ac_opf_solution_status(s, &status) == YBUS_OK && status == YBUS_OPF_CONVERGED);
    ybus_ac_opf_solution_free(s);
    ybus_case_free(case14);
}

int main(void) {
    ybus_case *three = three_bus();

    /* three_bus.m's generators have PQ capability curves, which the AC OPF can't use */
    ybus_ac_opf_solution *solution = NULL;
    CHECK(ybus_ac_opf(three, NULL, &solution) == YBUS_ERROR_UNSUPPORTED_GENERATOR);
    CHECK(solution == NULL && strstr(last_message(), "capability curve") != NULL);

    /* Without the curves, it solves */
    const double zero[2] = {0, 0};
    CHECK(ybus_case_set_gen_f64(three, YBUS_GEN_PC1_MW, zero, 2) == YBUS_OK);
    CHECK(ybus_case_set_gen_f64(three, YBUS_GEN_PC2_MW, zero, 2) == YBUS_OK);
    CHECK(ybus_ac_opf(three, NULL, &solution) == YBUS_OK);
    ybus_opf_status status = 9;
    CHECK(ybus_ac_opf_solution_status(solution, &status) == YBUS_OK);
    CHECK(status == YBUS_OPF_CONVERGED);
    double vm[3], vg[2], lam[3];
    CHECK(ybus_ac_opf_solution_get_bus(solution, YBUS_BUS_VM_PU, vm, 3) == YBUS_OK);
    CHECK(ybus_ac_opf_solution_get_gen(solution, YBUS_GEN_VG_PU, vg, 2) == YBUS_OK);
    CHECK(vg[0] == vm[0] && vg[1] == vm[1]);
    CHECK(ybus_ac_opf_solution_get_bus(solution, YBUS_BUS_LAM_P, lam, 3) == YBUS_OK);
    CHECK(lam[0] > 0);
    for (size_t i = 0; i < 3; i++) {
        CHECK(vm[i] >= 0.9 - 1e-6 && vm[i] <= 1.1 + 1e-6);
    }
    CHECK(ybus_ac_opf_solution_get_gen(solution, YBUS_GEN_PMAX_MW, vg, 2) ==
          YBUS_ERROR_INVALID_FIELD);
    ybus_ac_opf_solution_free(solution);
    ybus_ac_opf_solution_free(NULL);

    /* Options given, here the predictor-corrector */
    ybus_ac_opf_options *options = NULL;
    CHECK(ybus_ac_opf_options_new(YBUS_OPF_MEHROTRA, &options) == YBUS_OK);
    CHECK(ybus_ac_opf(three, options, &solution) == YBUS_OK);
    CHECK(ybus_ac_opf_solution_status(solution, &status) == YBUS_OK);
    CHECK(status == YBUS_OPF_CONVERGED);
    ybus_ac_opf_solution_free(solution);
    ybus_ac_opf_options_free(options);
    ybus_case_free(three);

    case14();
    return check_result();
}
