/* The last error message: a thread on which no call has failed has none, and the buffer is
 * left alone. The messages of real failures are tested from phase F3 on, once a call can fail
 * (loading a missing case file, for one). */

#include <string.h>

#include "check.h"
#include "ybus.h"

int main(void) {
    char buffer[16];
    memset(buffer, 'x', sizeof buffer);
    CHECK(ybus_last_error_message(NULL, 0) == 0);
    CHECK(ybus_last_error_message(buffer, sizeof buffer) == 0);
    CHECK(buffer[0] == 'x' && buffer[sizeof buffer - 1] == 'x');
    CHECK(YBUS_OK == 0);
    return check_result();
}
