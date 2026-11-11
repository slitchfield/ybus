/* The DC and AC optimal power flows' options: each solver's defaults, each setting set and
 * read back on either solver, and the values refused. */

#include <math.h>
#include <string.h>

#include "cases.h"
#include "check.h"

int main(void) {
    /* MATPOWER's defaults: MIPS, feastol 5e-6, LU; the AC OPF scales its cost by 1e-4 */
    ybus_dc_opf_options *dc = NULL;
    ybus_ac_opf_options *ac = NULL, *mehrotra = NULL;
    double value = 0;
    size_t iterations = 0;
    ybus_opf_solver solver = 0;
    ybus_kkt_factorization kkt = 0;
    CHECK(ybus_dc_opf_options_new(YBUS_OPF_MIPS, &dc) == YBUS_OK);
    CHECK(ybus_ac_opf_options_new(YBUS_OPF_MIPS, &ac) == YBUS_OK);
    CHECK(ybus_dc_opf_options_get_feastol(dc, &value) == YBUS_OK && value == 5e-6);
    CHECK(ybus_dc_opf_options_get_gradtol(dc, &value) == YBUS_OK && value == 1e-6);
    CHECK(ybus_dc_opf_options_get_cost_mult(dc, &value) == YBUS_OK && value == 1.0);
    CHECK(ybus_ac_opf_options_get_cost_mult(ac, &value) == YBUS_OK && value == 1e-4);
    CHECK(ybus_dc_opf_options_get_max_iterations(dc, &iterations) == YBUS_OK);
    CHECK(iterations == 150);
    CHECK(ybus_dc_opf_options_get_kkt_factorization(dc, &kkt) == YBUS_OK && kkt == YBUS_KKT_LU);
    CHECK(ybus_dc_opf_options_get_solver(dc, &solver) == YBUS_OK && solver == YBUS_OPF_MIPS);

    /* The predictor-corrector's tuned settings factor by LDLᵀ */
    CHECK(ybus_ac_opf_options_new(YBUS_OPF_MEHROTRA, &mehrotra) == YBUS_OK);
    CHECK(ybus_ac_opf_options_get_solver(mehrotra, &solver) == YBUS_OK);
    CHECK(solver == YBUS_OPF_MEHROTRA);
    CHECK(ybus_ac_opf_options_get_kkt_factorization(mehrotra, &kkt) == YBUS_OK);
    CHECK(kkt == YBUS_KKT_LDLT);
    CHECK(ybus_ac_opf_options_get_feastol(mehrotra, &value) == YBUS_OK && value == 5e-6);

    /* Each setting, set and read back, on either solver */
    CHECK(ybus_ac_opf_options_set_comptol(mehrotra, 1e-8) == YBUS_OK);
    CHECK(ybus_ac_opf_options_get_comptol(mehrotra, &value) == YBUS_OK && value == 1e-8);
    CHECK(ybus_dc_opf_options_set_costtol(dc, 2e-8) == YBUS_OK);
    CHECK(ybus_dc_opf_options_get_costtol(dc, &value) == YBUS_OK && value == 2e-8);
    CHECK(ybus_dc_opf_options_set_max_iterations(dc, 40) == YBUS_OK);
    CHECK(ybus_dc_opf_options_get_max_iterations(dc, &iterations) == YBUS_OK);
    CHECK(iterations == 40);
    CHECK(ybus_dc_opf_options_set_kkt_factorization(dc, YBUS_KKT_LDLT) == YBUS_OK);
    CHECK(ybus_dc_opf_options_get_kkt_factorization(dc, &kkt) == YBUS_OK);
    CHECK(kkt == YBUS_KKT_LDLT);

    /* Refused values change nothing */
    CHECK(ybus_dc_opf_options_set_feastol(dc, -1.0) == YBUS_ERROR_INVALID_VALUE);
    CHECK(strcmp(last_message(), "the feastol -1 isn't positive and finite") == 0);
    CHECK(ybus_ac_opf_options_set_cost_mult(ac, NAN) == YBUS_ERROR_INVALID_VALUE);
    CHECK(ybus_dc_opf_options_set_kkt_factorization(dc, 3) == YBUS_ERROR_INVALID_VALUE);
    CHECK(ybus_dc_opf_options_get_feastol(dc, &value) == YBUS_OK && value == 5e-6);
    CHECK(ybus_ac_opf_options_get_cost_mult(ac, &value) == YBUS_OK && value == 1e-4);
    ybus_dc_opf_options *out = dc;
    CHECK(ybus_dc_opf_options_new(0, &out) == YBUS_ERROR_INVALID_VALUE && out == NULL);
    CHECK(ybus_ac_opf_options_get_gradtol(NULL, &value) == YBUS_ERROR_NULL_ARGUMENT);

    ybus_dc_opf_options_free(dc);
    ybus_ac_opf_options_free(ac);
    ybus_ac_opf_options_free(mehrotra);
    ybus_dc_opf_options_free(NULL);
    ybus_ac_opf_options_free(NULL);
    return check_result();
}
