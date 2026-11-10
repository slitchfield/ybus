/* The C tests' cases: three_bus.m, the case the Rust tests also use, loaded from the
 * directory that run.sh names in YBUS_TEST_DATA, and the last error message. Every function is
 * static inline, so a test that doesn't use one isn't warned about it. */

#ifndef YBUS_TEST_CASES_H
#define YBUS_TEST_CASES_H

#include <stdio.h>
#include <stdlib.h>

#include "ybus.h"

/* The path of a file in the test data directory, in a buffer reused by the next call */
static inline const char *test_data(const char *name) {
    static char path[4096];
    const char *dir = getenv("YBUS_TEST_DATA");
    if (dir == NULL) {
        fprintf(stderr, "YBUS_TEST_DATA isn't set: run the tests with tests/c/run.sh\n");
        exit(2);
    }
    snprintf(path, sizeof path, "%s/%s", dir, name);
    return path;
}

/* The thread's last error message, or "" if there's none, in a buffer reused by the next call */
static inline const char *last_message(void) {
    static char message[1024];
    size_t needed = ybus_last_error_message(message, sizeof message);
    if (needed == 0 || needed > sizeof message) {
        message[0] = '\0';
    }
    return message;
}

/* The three-bus case, which the test frees; exits the test if it can't be loaded */
static inline ybus_case *three_bus(void) {
    ybus_case *loaded = NULL;
    if (ybus_case_from_path(test_data("three_bus.m"), &loaded) != YBUS_OK) {
        fprintf(stderr, "three_bus.m didn't load: %s\n", last_message());
        exit(2);
    }
    return loaded;
}

/* A whole file, in a buffer the caller frees, with its length in *len; exits if it can't be
 * read */
static inline unsigned char *read_file(const char *path, size_t *len) {
    FILE *file = fopen(path, "rb");
    unsigned char *bytes = NULL;
    long size;
    if (file == NULL || fseek(file, 0, SEEK_END) != 0 || (size = ftell(file)) < 0 ||
        fseek(file, 0, SEEK_SET) != 0 || (bytes = malloc((size_t)size + 1)) == NULL ||
        fread(bytes, 1, (size_t)size, file) != (size_t)size) {
        fprintf(stderr, "%s couldn't be read\n", path);
        exit(2);
    }
    fclose(file);
    *len = (size_t)size;
    return bytes;
}

#endif /* YBUS_TEST_CASES_H */
