/* The C tests' one check macro, so they need no framework: CHECK reports a failed condition
 * with its file and line and counts it, and main returns check_result(), which is non-zero if
 * any check failed. */

#ifndef YBUS_TEST_CHECK_H
#define YBUS_TEST_CHECK_H

#include <stdio.h>

static int check_failures = 0;

#define CHECK(condition)                                                                   \
    do {                                                                                   \
        if (!(condition)) {                                                                \
            fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
            check_failures++;                                                              \
        }                                                                                  \
    } while (0)

static int check_result(void) {
    if (check_failures > 0) {
        fprintf(stderr, "%d check(s) failed\n", check_failures);
        return 1;
    }
    return 0;
}

#endif /* YBUS_TEST_CHECK_H */
