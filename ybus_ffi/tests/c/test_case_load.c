/* Loading cases from a file and from bytes, copying and freeing them, and the errors C sees
 * when a case can't be loaded. */

#include <stdlib.h>
#include <string.h>

#include "cases.h"
#include "check.h"

int main(void) {
    /* From a file, from its bytes, and copied; then freed, with NULL freed too */
    ybus_case *from_file = NULL, *from_bytes = NULL, *copy = NULL;
    size_t len = 0;
    unsigned char *bytes = read_file(test_data("three_bus.m"), &len);
    CHECK(ybus_case_from_path(test_data("three_bus.m"), &from_file) == YBUS_OK);
    CHECK(ybus_case_from_bytes(bytes, len, &from_bytes) == YBUS_OK);
    CHECK(ybus_case_clone(from_file, &copy) == YBUS_OK);
    CHECK(from_file != NULL && from_bytes != NULL && copy != NULL);
    free(bytes);
    ybus_case_free(from_file);
    ybus_case_free(from_bytes);
    ybus_case_free(copy);
    ybus_case_free(NULL);

    /* Each error leaves *out NULL, even when it held a case before */
    ybus_case *held = three_bus();
    ybus_case *out = held;
    CHECK(ybus_case_from_path("/nonexistent/ybus/case.m", &out) == YBUS_ERROR_NO_CASE_FILE);
    CHECK(out == NULL);
    CHECK(strstr(last_message(), "/nonexistent/ybus/case.m") != NULL);

    const char malformed[] = "mpc.version = '2';\nmpc.bus = [1 2;\n";
    out = held;
    CHECK(ybus_case_from_bytes((const uint8_t *)malformed, strlen(malformed), &out) ==
          YBUS_ERROR_PARSE);
    CHECK(out == NULL);

    CHECK(ybus_case_from_path("/tmp/\xff.m", &out) == YBUS_ERROR_INVALID_UTF8);
    CHECK(strcmp(last_message(), "`path` isn't valid UTF-8") == 0);

    CHECK(ybus_case_from_path(test_data("three_bus.m"), NULL) == YBUS_ERROR_NULL_ARGUMENT);
    CHECK(strcmp(last_message(), "`out` is null") == 0);
    CHECK(ybus_case_clone(NULL, &out) == YBUS_ERROR_NULL_ARGUMENT);
    CHECK(strcmp(last_message(), "`case` is null") == 0);
    ybus_case_free(held);

    return check_result();
}
