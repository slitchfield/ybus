/* A case's columns, read and written whole by field code, and its bus names, read one at a
 * time. */

#include <stdint.h>
#include <string.h>

#include "cases.h"
#include "check.h"

int main(void) {
    ybus_case *three = three_bus();

    /* A double column and the int64_t ones, against three_bus.m's text */
    double pd[3] = {0};
    CHECK(ybus_case_get_bus_f64(three, YBUS_BUS_PD_MW, pd, 3) == YBUS_OK);
    CHECK(pd[0] == 0.0 && pd[1] == 50.0 && pd[2] == 80.0);
    int64_t types[3] = {0};
    CHECK(ybus_case_get_bus_i64(three, YBUS_BUS_TYPE, types, 3) == YBUS_OK);
    CHECK(types[0] == 3 && types[1] == 2 && types[2] == 1);
    double apf[2] = {0};
    CHECK(ybus_case_get_gen_f64(three, YBUS_GEN_APF, apf, 2) == YBUS_OK);
    CHECK(apf[0] == 21.0 && apf[1] == 2.1);
    int64_t status[3] = {0};
    CHECK(ybus_case_get_branch_i64(three, YBUS_BRANCH_STATUS, status, 3) == YBUS_OK);
    CHECK(status[0] == 1 && status[1] == 1 && status[2] == 0);
    int64_t to[3] = {0};
    CHECK(ybus_case_get_branch_i64(three, YBUS_BRANCH_TO_LABEL, to, 3) == YBUS_OK);
    CHECK(to[0] == 2 && to[1] == 3 && to[2] == 3);
    double x[3] = {0};
    CHECK(ybus_case_get_branch_f64(three, YBUS_BRANCH_X, x, 3) == YBUS_OK);
    CHECK(x[0] == 0.1 && x[1] == 0.2 && x[2] == 0.25);

    /* A code the table doesn't have, or of the other type */
    CHECK(ybus_case_get_bus_f64(three, 14, pd, 3) == YBUS_ERROR_INVALID_FIELD);
    CHECK(ybus_case_get_bus_f64(three, YBUS_BUS_TYPE, pd, 3) == YBUS_ERROR_INVALID_FIELD);
    CHECK(strcmp(last_message(), "2 isn't a double field of the bus table") == 0);
    CHECK(ybus_case_get_gen_i64(three, YBUS_GEN_PG_MW, types, 3) == YBUS_ERROR_INVALID_FIELD);

    /* Too short a buffer is refused untouched; a longer one keeps its tail */
    double qd[5] = {-1, -1, -1, -1, -1};
    CHECK(ybus_case_get_bus_f64(three, YBUS_BUS_QD_MVAR, qd, 2) == YBUS_ERROR_BUFFER_TOO_SMALL);
    CHECK(qd[0] == -1 && qd[1] == -1);
    CHECK(ybus_case_get_bus_f64(three, YBUS_BUS_QD_MVAR, qd, 5) == YBUS_OK);
    CHECK(qd[0] == 0 && qd[1] == 10 && qd[2] == 20 && qd[3] == -1 && qd[4] == -1);

    /* Bus names: the size first, then the name */
    const char *names[3] = {"Alpha", "Beta", "Gamma"};
    for (size_t row = 0; row < 3; row++) {
        size_t needed = 0;
        char name[16];
        CHECK(ybus_case_bus_name(three, row, NULL, 0, &needed) == YBUS_OK);
        CHECK(needed == strlen(names[row]) + 1);
        CHECK(ybus_case_bus_name(three, row, name, needed, NULL) == YBUS_OK);
        CHECK(strcmp(name, names[row]) == 0);
    }
    char short_name[3] = {'x', 'x', 'x'};
    size_t needed = 0;
    CHECK(ybus_case_bus_name(three, 0, short_name, sizeof short_name, &needed) ==
          YBUS_ERROR_BUFFER_TOO_SMALL);
    CHECK(needed == 6 && short_name[0] == 'x');
    char name[16];
    CHECK(ybus_case_bus_name(three, 3, name, sizeof name, NULL) == YBUS_ERROR_ROW_OUT_OF_RANGE);

    /* Set a load and a status column, and read them back */
    const double loads[3] = {10.0, 20.0, 30.0};
    CHECK(ybus_case_set_bus_f64(three, YBUS_BUS_PD_MW, loads, 3) == YBUS_OK);
    CHECK(ybus_case_get_bus_f64(three, YBUS_BUS_PD_MW, pd, 3) == YBUS_OK);
    CHECK(pd[0] == 10.0 && pd[1] == 20.0 && pd[2] == 30.0);
    const int64_t switched[3] = {0, 1, 1};
    CHECK(ybus_case_set_branch_i64(three, YBUS_BRANCH_STATUS, switched, 3) == YBUS_OK);
    CHECK(ybus_case_get_branch_i64(three, YBUS_BRANCH_STATUS, status, 3) == YBUS_OK);
    CHECK(status[0] == 0 && status[1] == 1 && status[2] == 1);

    /* A value the field can't hold leaves the whole column as it was */
    const int64_t bad_types[3] = {1, 2, 7};
    CHECK(ybus_case_set_bus_i64(three, YBUS_BUS_TYPE, bad_types, 3) == YBUS_ERROR_INVALID_VALUE);
    CHECK(strcmp(last_message(), "bus row 2's value 7 for field 2 isn't a bus type, 1 to 4") == 0);
    CHECK(ybus_case_get_bus_i64(three, YBUS_BUS_TYPE, types, 3) == YBUS_OK);
    CHECK(types[0] == 3 && types[1] == 2 && types[2] == 1);

    /* One entry per row, no more and no fewer */
    CHECK(ybus_case_set_gen_f64(three, YBUS_GEN_PG_MW, loads, 3) == YBUS_ERROR_LENGTH_MISMATCH);
    CHECK(ybus_case_set_gen_f64(three, YBUS_GEN_PG_MW, loads, 2) == YBUS_OK);

    ybus_case_free(three);
    return check_result();
}
