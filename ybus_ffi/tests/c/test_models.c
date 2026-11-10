/* DC and AC models: built once and solved for several load scenarios, each the same, bit for
 * bit, as the power flow of the same case; an AC re-stamp after a reactance changes; and a
 * case of another topology refused. */

#include <string.h>

#include "cases.h"
#include "check.h"

/* Every column of two DC solutions is the same, bit for bit */
static int dc_same(const ybus_dc_solution *a, const ybus_dc_solution *b) {
    double x[3], y[3];
    return ybus_dc_solution_get_bus(a, YBUS_BUS_VA_DEG, x, 3) == YBUS_OK &&
           ybus_dc_solution_get_bus(b, YBUS_BUS_VA_DEG, y, 3) == YBUS_OK &&
           memcmp(x, y, sizeof x) == 0 &&
           ybus_dc_solution_get_gen(a, YBUS_GEN_PG_MW, x, 2) == YBUS_OK &&
           ybus_dc_solution_get_gen(b, YBUS_GEN_PG_MW, y, 2) == YBUS_OK &&
           memcmp(x, y, 2 * sizeof x[0]) == 0 &&
           ybus_dc_solution_get_branch(a, YBUS_BRANCH_PF_MW, x, 3) == YBUS_OK &&
           ybus_dc_solution_get_branch(b, YBUS_BRANCH_PF_MW, y, 3) == YBUS_OK &&
           memcmp(x, y, sizeof x) == 0;
}

/* Every column of two AC solutions is the same, bit for bit */
static int ac_same(const ybus_ac_solution *a, const ybus_ac_solution *b) {
    const ybus_field buses[] = {YBUS_BUS_VM_PU, YBUS_BUS_VA_DEG};
    const ybus_field gens[] = {YBUS_GEN_PG_MW, YBUS_GEN_QG_MVAR};
    const ybus_field branches[] = {YBUS_BRANCH_PF_MW, YBUS_BRANCH_QF_MVAR, YBUS_BRANCH_PT_MW,
                                   YBUS_BRANCH_QT_MVAR};
    double x[3], y[3];
    for (size_t i = 0; i < 2; i++) {
        if (ybus_ac_solution_get_bus(a, buses[i], x, 3) != YBUS_OK ||
            ybus_ac_solution_get_bus(b, buses[i], y, 3) != YBUS_OK ||
            memcmp(x, y, sizeof x) != 0 ||
            ybus_ac_solution_get_gen(a, gens[i], x, 2) != YBUS_OK ||
            ybus_ac_solution_get_gen(b, gens[i], y, 2) != YBUS_OK ||
            memcmp(x, y, 2 * sizeof x[0]) != 0) {
            return 0;
        }
    }
    for (size_t i = 0; i < 4; i++) {
        if (ybus_ac_solution_get_branch(a, branches[i], x, 3) != YBUS_OK ||
            ybus_ac_solution_get_branch(b, branches[i], y, 3) != YBUS_OK ||
            memcmp(x, y, sizeof x) != 0) {
            return 0;
        }
    }
    return 1;
}

/* Scales every bus's load in the case */
static void scale_loads(ybus_case *scenario, double scale) {
    double pd[3], qd[3];
    CHECK(ybus_case_get_bus_f64(scenario, YBUS_BUS_PD_MW, pd, 3) == YBUS_OK);
    CHECK(ybus_case_get_bus_f64(scenario, YBUS_BUS_QD_MVAR, qd, 3) == YBUS_OK);
    for (size_t i = 0; i < 3; i++) {
        pd[i] *= scale;
        qd[i] *= scale;
    }
    CHECK(ybus_case_set_bus_f64(scenario, YBUS_BUS_PD_MW, pd, 3) == YBUS_OK);
    CHECK(ybus_case_set_bus_f64(scenario, YBUS_BUS_QD_MVAR, qd, 3) == YBUS_OK);
}

int main(void) {
    ybus_case *scenario = three_bus();
    ybus_dc_model *dc = NULL;
    ybus_ac_model *ac = NULL;
    CHECK(ybus_dc_model_new(scenario, &dc) == YBUS_OK);
    CHECK(ybus_ac_model_new(scenario, &ac) == YBUS_OK);

    /* Three load scenarios on the same models */
    const double scales[3] = {1.0, 0.75, 1.25};
    for (size_t i = 0; i < 3; i++) {
        scale_loads(scenario, scales[i]);
        ybus_dc_solution *on_model = NULL, *direct = NULL;
        CHECK(ybus_dc_model_solve(dc, scenario, NULL, &on_model) == YBUS_OK);
        CHECK(ybus_dc_power_flow(scenario, NULL, &direct) == YBUS_OK);
        CHECK(dc_same(on_model, direct));
        ybus_dc_solution_free(on_model);
        ybus_dc_solution_free(direct);

        ybus_ac_solution *ac_on_model = NULL, *ac_direct = NULL;
        CHECK(ybus_ac_model_solve(ac, scenario, NULL, &ac_on_model) == YBUS_OK);
        CHECK(ybus_ac_power_flow(scenario, NULL, &ac_direct) == YBUS_OK);
        CHECK(ac_same(ac_on_model, ac_direct));
        ybus_ac_solution_free(ac_on_model);
        ybus_ac_solution_free(ac_direct);
    }

    /* A changed reactance, re-stamped into the AC model */
    double x[3];
    CHECK(ybus_case_get_branch_f64(scenario, YBUS_BRANCH_X, x, 3) == YBUS_OK);
    x[1] = 0.3;
    CHECK(ybus_case_set_branch_f64(scenario, YBUS_BRANCH_X, x, 3) == YBUS_OK);
    CHECK(ybus_ac_model_restamp(ac, scenario) == YBUS_OK);
    ybus_ac_solution *restamped = NULL, *direct = NULL;
    CHECK(ybus_ac_model_solve(ac, scenario, NULL, &restamped) == YBUS_OK);
    CHECK(ybus_ac_power_flow(scenario, NULL, &direct) == YBUS_OK);
    CHECK(ac_same(restamped, direct));
    ybus_ac_solution_free(restamped);
    ybus_ac_solution_free(direct);

    /* A branch switched out changes the topology, which needs new models */
    const int64_t switched[3] = {0, 1, 1};
    CHECK(ybus_case_set_branch_i64(scenario, YBUS_BRANCH_STATUS, switched, 3) == YBUS_OK);
    ybus_dc_solution *refused = NULL;
    CHECK(ybus_dc_model_solve(dc, scenario, NULL, &refused) == YBUS_ERROR_TOPOLOGY_MISMATCH);
    CHECK(refused == NULL);
    CHECK(strstr(last_message(), "branch row 1 has a different status") != NULL);
    CHECK(ybus_ac_model_restamp(ac, scenario) == YBUS_ERROR_TOPOLOGY_MISMATCH);
    CHECK(ybus_dc_model_solve(NULL, scenario, NULL, &refused) == YBUS_ERROR_NULL_ARGUMENT);
    CHECK(strcmp(last_message(), "`model` is null") == 0);

    ybus_dc_model_free(dc);
    ybus_ac_model_free(ac);
    ybus_dc_model_free(NULL);
    ybus_ac_model_free(NULL);
    ybus_case_free(scenario);
    return check_result();
}
