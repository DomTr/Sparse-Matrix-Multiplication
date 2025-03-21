#include <stdint.h>
#include <stdio.h>
#include "matr_mult_ellpack_v1.h"
/*
 * function to multiply two EllpackMatrix and save it to the result pointer
 */
void matr_mult_ellpack_v1(const void *a, const void *b, void *result)
{
    sequential_multiplication(a, b, result);
}
