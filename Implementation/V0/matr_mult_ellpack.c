#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "matr_mult_ellpack.h"
#include "../utils.h"

/*
 * Main method to multiply two EllpackMatrix and save it to the result pointer
 */
void matr_mult_ellpack(const void *a, const void *b, void *result)
{
    EllpackMatrix *a_matrix = (EllpackMatrix *)a;
    EllpackMatrix *b_matrix = (EllpackMatrix *)b;
    float **result_matrix = (float **)result;

    if (a_matrix->cols != b_matrix->rows)
    {
        fprintf(stderr, ERR_INVALID_MATRIX_DIMENSIONS, a_matrix->cols, b_matrix->rows);
        free_matrix_array(a_matrix->rows, result);
        free_ellpack_matrix((EllpackMatrix *)a);
        free_ellpack_matrix((EllpackMatrix *)b);
        exit(EXIT_FAILURE);
    }

    for (uint64_t a_row = 0; a_row < a_matrix->rows; a_row++)
    {
        for (uint64_t a_ellpack_col = 0; a_ellpack_col < a_matrix->ellpack_cols; a_ellpack_col++)
        {
            float a_val = a_matrix->values[a_row][a_ellpack_col];
            // break the row iteration once we encountered a zero entry
            if (a_val == 0)
            {
                break;
            }
            uint64_t a_col = a_matrix->indices[a_row][a_ellpack_col];

            // b_row_vector contains values in b_matrix that a_val can multiply
            // b_row_indices contains the column index of the values in b_matrix
            float *b_row_vector = b_matrix->values[a_col];
            uint64_t *b_row_indices = b_matrix->indices[a_col];

            // result row would be a_row
            // result col would be b_col
            scalar_multiplication_v0(a_val, b_row_vector, b_row_indices, result_matrix[a_row], b_matrix->ellpack_cols);
        }
    }
}

/*
 * scalar multiplication a * b_row_vector
 * multiplication result is accumulated to each entry in result_row_vector
 */
void scalar_multiplication_v0(float a_val, const float *b_row_vector, const uint64_t *b_row_indices, float *result_row_vector, uint64_t b_ellpack_cols)
{
    for (uint64_t i = 0; i < b_ellpack_cols; i++)
    {
        float b_val = b_row_vector[i];
        if (b_val == 0)
        {
            break;
        }
        uint64_t b_col = b_row_indices[i];
        result_row_vector[b_col] += a_val * b_val;
    }
}
