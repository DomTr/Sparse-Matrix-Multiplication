#ifndef FINAL_MATR_MULT_ELLPACK_V2_H
#define FINAL_MATR_MULT_ELLPACK_V2_H

#include <pthread.h>
#include "../utils.h"
#include "../optimizations.h"

//ThreadData struct
typedef struct{
    uint64_t start;
    uint64_t end;
    EllpackMatrix *a_matrix;
    EllpackMatrix *b_matrix;
    float **result;
} ThreadData;

void *matr_mult_one_thread(void *arg);

void parallel_multiplication(const void *a, const void *b, void *result);

void matr_mult_ellpack_v2(const void *a, const void *b, void *result);

#endif
