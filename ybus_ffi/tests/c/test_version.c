/* The library's version: the run-time library's is the header's, and the header's encoding
 * agrees with its parts. */

#include "check.h"
#include "ybus.h"

int main(void) {
    CHECK(ybus_version() == YBUS_VERSION);
    CHECK(YBUS_VERSION ==
          YBUS_VERSION_MAJOR * 1000000 + YBUS_VERSION_MINOR * 1000 + YBUS_VERSION_PATCH);
    /* The encoding holds only while minor and patch stay below 1000 */
    CHECK(YBUS_VERSION_MINOR < 1000 && YBUS_VERSION_PATCH < 1000);
    return check_result();
}
