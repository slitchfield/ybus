/* The C tests' cases: three_bus.m, the case the Rust tests also use, loaded from the
 * directory that run.sh names in YBUS_TEST_DATA; MATPOWER's cases and their committed
 * references, from YBUS_MATPOWER_DATA and YBUS_REFERENCE_DATA; and the last error message.
 * Every function is static inline, so a test that doesn't use one isn't warned about it. */

#ifndef YBUS_TEST_CASES_H
#define YBUS_TEST_CASES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

/* The path of MATPOWER's <name>.m, in a buffer reused by the next call, or NULL if it isn't
 * fetched: the test skips it, saying so, unless YBUS_REQUIRE_CASES=1, when it fails instead */
static inline const char *matpower_case(const char *name) {
    static char path[4096];
    const char *dir = getenv("YBUS_MATPOWER_DATA");
    const char *require = getenv("YBUS_REQUIRE_CASES");
    FILE *file = NULL;
    if (dir != NULL) {
        snprintf(path, sizeof path, "%s/%s.m", dir, name);
        file = fopen(path, "rb");
    }
    if (file != NULL) {
        fclose(file);
        return path;
    }
    if (require != NULL && strcmp(require, "1") == 0) {
        fprintf(stderr, "MATPOWER's %s.m is missing: run scripts/fetch-matpower.sh\n", name);
        exit(1);
    }
    fprintf(stderr, "skipped: MATPOWER's %s.m isn't fetched\n", name);
    return NULL;
}

/* Reads `column` of the committed reference <case>_<run>_<table>.csv into out, at most max
 * values, and returns how many it read; exits if the file or the column isn't there */
static inline size_t read_reference(const char *case_name, const char *run, const char *table,
                                    const char *column, double *out, size_t max) {
    char path[4096], line[4096];
    const char *dir = getenv("YBUS_REFERENCE_DATA");
    snprintf(path, sizeof path, "%s/%s_%s_%s.csv", dir == NULL ? "." : dir, case_name, run, table);
    FILE *file = fopen(path, "r");
    if (file == NULL || fgets(line, sizeof line, file) == NULL) {
        fprintf(stderr, "%s couldn't be read\n", path);
        exit(2);
    }
    /* The column's position in the header */
    int index = -1, position = 0;
    for (char *name = strtok(line, ",\r\n"); name != NULL; name = strtok(NULL, ",\r\n")) {
        if (strcmp(name, column) == 0) {
            index = position;
        }
        position++;
    }
    if (index < 0) {
        fprintf(stderr, "%s has no column %s\n", path, column);
        exit(2);
    }
    size_t count = 0;
    while (count < max && fgets(line, sizeof line, file) != NULL) {
        char *field = strtok(line, ",\r\n");
        for (int i = 0; i < index && field != NULL; i++) {
            field = strtok(NULL, ",\r\n");
        }
        if (field == NULL) {
            fprintf(stderr, "%s has a short row\n", path);
            exit(2);
        }
        out[count++] = strtod(field, NULL);
    }
    fclose(file);
    return count;
}

/* Whether every one of n values is within tolerance of its reference */
static inline int all_close(const double *values, const double *reference, size_t n,
                            double tolerance) {
    for (size_t i = 0; i < n; i++) {
        double difference = values[i] - reference[i];
        if (difference > tolerance || difference < -tolerance) {
            fprintf(stderr, "entry %zu: %.17g against %.17g\n", i, values[i], reference[i]);
            return 0;
        }
    }
    return 1;
}

#endif /* YBUS_TEST_CASES_H */
