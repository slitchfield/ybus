/* A case's sizes and scalars: its table lengths, MVA base, and whether it has costs and bus
 * names. */

#include <stdbool.h>
#include <string.h>

#include "cases.h"
#include "check.h"

int main(void) {
    ybus_case *three = three_bus();
    size_t buses = 0, generators = 0, branches = 0;
    double base = 0.0;
    bool costs = false, names = false;
    CHECK(ybus_case_bus_count(three, &buses) == YBUS_OK && buses == 3);
    CHECK(ybus_case_gen_count(three, &generators) == YBUS_OK && generators == 2);
    CHECK(ybus_case_branch_count(three, &branches) == YBUS_OK && branches == 3);
    CHECK(ybus_case_base_mva(three, &base) == YBUS_OK && base == 100.0);
    CHECK(ybus_case_has_costs(three, &costs) == YBUS_OK && costs);
    CHECK(ybus_case_has_bus_names(three, &names) == YBUS_OK && names);

    CHECK(ybus_case_set_base_mva(three, 50.0) == YBUS_OK);
    CHECK(ybus_case_base_mva(three, &base) == YBUS_OK && base == 50.0);

    CHECK(ybus_case_bus_count(NULL, &buses) == YBUS_ERROR_NULL_ARGUMENT);
    CHECK(strcmp(last_message(), "`case` is null") == 0);
    CHECK(ybus_case_base_mva(three, NULL) == YBUS_ERROR_NULL_ARGUMENT);
    CHECK(strcmp(last_message(), "`out` is null") == 0);
    CHECK(ybus_case_set_base_mva(NULL, 1.0) == YBUS_ERROR_NULL_ARGUMENT);

    ybus_case_free(three);
    return check_result();
}
