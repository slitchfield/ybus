/* The DC and AC power flow options: their defaults, each setting set and read back, and the
 * values refused. */

#include <math.h>
#include <string.h>

#include "cases.h"
#include "check.h"

int main(void) {
    /* DC: AMD by default, and each ordering by its code */
    ybus_dc_options *dc = NULL;
    ybus_ordering ordering = 0;
    CHECK(ybus_dc_options_new(&dc) == YBUS_OK);
    CHECK(ybus_dc_options_get_ordering(dc, &ordering) == YBUS_OK);
    CHECK(ordering == YBUS_ORDERING_AMD);
    CHECK(ybus_dc_options_set_ordering(dc, YBUS_ORDERING_RCM) == YBUS_OK);
    CHECK(ybus_dc_options_get_ordering(dc, &ordering) == YBUS_OK);
    CHECK(ordering == YBUS_ORDERING_RCM);
    CHECK(ybus_dc_options_set_ordering(dc, 9) == YBUS_ERROR_INVALID_VALUE);
    CHECK(strcmp(last_message(), "9 isn't an ordering, 1 to 4") == 0);
    ybus_dc_options_free(dc);
    ybus_dc_options_free(NULL);

    /* AC: MATPOWER's defaults for each method */
    ybus_ac_options *newton = NULL, *fast = NULL;
    double tolerance = 0;
    size_t iterations = 0;
    ybus_ac_method method = 0;
    ybus_reactive_limits limits = -1;
    CHECK(ybus_ac_options_new(YBUS_AC_NEWTON, &newton) == YBUS_OK);
    CHECK(ybus_ac_options_new(YBUS_AC_FAST_DECOUPLED_XB, &fast) == YBUS_OK);
    CHECK(ybus_ac_options_get_tolerance(newton, &tolerance) == YBUS_OK && tolerance == 1e-8);
    CHECK(ybus_ac_options_get_max_iterations(newton, &iterations) == YBUS_OK && iterations == 10);
    CHECK(ybus_ac_options_get_max_iterations(fast, &iterations) == YBUS_OK && iterations == 30);
    CHECK(ybus_ac_options_get_method(fast, &method) == YBUS_OK);
    CHECK(method == YBUS_AC_FAST_DECOUPLED_XB);
    CHECK(ybus_ac_options_get_reactive_limits(newton, &limits) == YBUS_OK);
    CHECK(limits == YBUS_REACTIVE_LIMITS_IGNORE);

    /* Each setting, set and read back; a method change keeps the iteration limit */
    CHECK(ybus_ac_options_set_method(newton, YBUS_AC_FAST_DECOUPLED_BX) == YBUS_OK);
    CHECK(ybus_ac_options_set_tolerance(newton, 1e-6) == YBUS_OK);
    CHECK(ybus_ac_options_set_reactive_limits(newton, YBUS_REACTIVE_LIMITS_ALL_AT_ONCE) ==
          YBUS_OK);
    CHECK(ybus_ac_options_get_method(newton, &method) == YBUS_OK);
    CHECK(method == YBUS_AC_FAST_DECOUPLED_BX);
    CHECK(ybus_ac_options_get_max_iterations(newton, &iterations) == YBUS_OK && iterations == 10);
    CHECK(ybus_ac_options_set_max_iterations(newton, 40) == YBUS_OK);
    CHECK(ybus_ac_options_get_max_iterations(newton, &iterations) == YBUS_OK && iterations == 40);
    CHECK(ybus_ac_options_get_tolerance(newton, &tolerance) == YBUS_OK && tolerance == 1e-6);
    CHECK(ybus_ac_options_get_reactive_limits(newton, &limits) == YBUS_OK);
    CHECK(limits == YBUS_REACTIVE_LIMITS_ALL_AT_ONCE);

    /* Refused values change nothing */
    CHECK(ybus_ac_options_set_tolerance(newton, 0.0) == YBUS_ERROR_INVALID_VALUE);
    CHECK(ybus_ac_options_set_tolerance(newton, NAN) == YBUS_ERROR_INVALID_VALUE);
    CHECK(ybus_ac_options_set_method(newton, 0) == YBUS_ERROR_INVALID_VALUE);
    CHECK(ybus_ac_options_set_reactive_limits(newton, 3) == YBUS_ERROR_INVALID_VALUE);
    CHECK(ybus_ac_options_get_tolerance(newton, &tolerance) == YBUS_OK && tolerance == 1e-6);
    ybus_ac_options *out = newton;
    CHECK(ybus_ac_options_new(7, &out) == YBUS_ERROR_INVALID_VALUE && out == NULL);
    CHECK(ybus_ac_options_get_tolerance(NULL, &tolerance) == YBUS_ERROR_NULL_ARGUMENT);
    CHECK(strcmp(last_message(), "`options` is null") == 0);

    ybus_ac_options_free(newton);
    ybus_ac_options_free(fast);
    ybus_ac_options_free(NULL);
    return check_result();
}
