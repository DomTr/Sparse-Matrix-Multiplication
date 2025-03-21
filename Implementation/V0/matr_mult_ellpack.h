#ifndef FINAL_MATR_MULT_ELLPACK_H
#define FINAL_MATR_MULT_ELLPACK_H

void matr_mult_ellpack(const void *a, const void *b, void *result);

void scalar_multiplication_v0(float a_val, const float *b_row_vector, const uint64_t *b_row_indices, float *result_row_vector, uint64_t b_ellpack_cols);

#endif
