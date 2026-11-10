/* A case's generator costs, read and set whole as flat arrays with offsets. */

#include <string.h>

#include "cases.h"
#include "check.h"

int main(void) {
    ybus_case *three = three_bus();

    /* three_bus.m's two rows: a quadratic, and a piecewise linear cost through three points */
    size_t rows = 0, count = 0;
    CHECK(ybus_case_cost_sizes(three, &rows, &count) == YBUS_OK);
    CHECK(rows == 2 && count == 9);
    int32_t model[2];
    double startup[2], shutdown[2], values[9];
    size_t offsets[3];
    CHECK(ybus_case_get_costs(three, model, startup, shutdown, offsets, 2, values, 9) == YBUS_OK);
    CHECK(model[0] == YBUS_COST_POLYNOMIAL && model[1] == YBUS_COST_PIECEWISE_LINEAR);
    CHECK(startup[0] == 1500 && shutdown[0] == 100 && startup[1] == 0 && shutdown[1] == 0);
    CHECK(offsets[0] == 0 && offsets[1] == 3 && offsets[2] == 9);
    CHECK(values[0] == 0.11 && values[1] == 5 && values[2] == 150);
    CHECK(values[3] == 0 && values[4] == 0 && values[7] == 80 && values[8] == 4400);

    /* New costs, both linear, read back as set */
    const int32_t new_model[2] = {YBUS_COST_POLYNOMIAL, YBUS_COST_PIECEWISE_LINEAR};
    const double new_startup[2] = {0, 0}, new_shutdown[2] = {0, 0};
    const size_t new_offsets[3] = {0, 2, 6};
    const double new_values[6] = {20, 0, 0, 0, 100, 2500};
    CHECK(ybus_case_set_costs(three, new_model, new_startup, new_shutdown, new_offsets, 2,
                              new_values, 6) == YBUS_OK);
    CHECK(ybus_case_cost_sizes(three, &rows, &count) == YBUS_OK && rows == 2 && count == 6);
    CHECK(ybus_case_get_costs(three, model, startup, shutdown, offsets, 2, values, 9) == YBUS_OK);
    CHECK(offsets[1] == 2 && offsets[2] == 6 && values[0] == 20 && values[5] == 2500);

    /* A malformed table changes nothing */
    const int32_t bad_model[2] = {YBUS_COST_POLYNOMIAL, 3};
    CHECK(ybus_case_set_costs(three, bad_model, new_startup, new_shutdown, new_offsets, 2,
                              new_values, 6) == YBUS_ERROR_INVALID_VALUE);
    const size_t odd[3] = {0, 1, 6};
    CHECK(ybus_case_set_costs(three, new_model, new_startup, new_shutdown, odd, 2, new_values,
                              6) == YBUS_ERROR_INVALID_VALUE);
    const size_t short_offsets[3] = {0, 2, 5};
    CHECK(ybus_case_set_costs(three, new_model, new_startup, new_shutdown, short_offsets, 2,
                              new_values, 6) == YBUS_ERROR_LENGTH_MISMATCH);
    CHECK(ybus_case_get_costs(three, model, startup, shutdown, offsets, 2, values, 9) == YBUS_OK);
    CHECK(model[1] == YBUS_COST_PIECEWISE_LINEAR && offsets[2] == 6 && values[5] == 2500);

    /* Too little room, and then no costs at all */
    CHECK(ybus_case_get_costs(three, model, startup, shutdown, offsets, 2, values, 5) ==
          YBUS_ERROR_BUFFER_TOO_SMALL);
    CHECK(ybus_case_set_costs(three, NULL, NULL, NULL, NULL, 0, NULL, 0) == YBUS_OK);
    bool costs = true;
    CHECK(ybus_case_has_costs(three, &costs) == YBUS_OK && !costs);
    CHECK(ybus_case_cost_sizes(three, &rows, &count) == YBUS_OK && rows == 0 && count == 0);
    CHECK(ybus_case_get_costs(three, model, startup, shutdown, offsets, 2, values, 9) ==
          YBUS_ERROR_MISSING_DATA);
    CHECK(strcmp(last_message(), "the case has no costs") == 0);

    ybus_case_free(three);
    return check_result();
}
