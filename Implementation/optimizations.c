#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "immintrin.h" // for simd
#include "optimizations.h"

/*
 * simd instructions are used to speed up multiplication itself
 * the last 8 instructions in the loop add sums to the result_row_vector. This is an unrolled loop
 * _mm_set_ps1(a_val) - all floats in m128 are set to value a_val
 * _mm_load_ps(address) - 4 floats are loaded into the register
 */
void scalar_multiplication_simd(float a_val, float *b_row_vector, const uint64_t *b_row_indices, float *result_row_vector, uint64_t b_ellpack_cols)
{
    size_t i = 0;
    __m128 sum1;
    __m128 sum2;
    for (; i < (b_ellpack_cols & ~7ul); i += 8)
    {
        sum1 = _mm_set_ps1(a_val); // set 4 values to a
        sum2 = _mm_set_ps1(a_val);
        sum1 *= _mm_load_ps(&b_row_vector[i]);    // a * b[i], a * b[i+1], a * b[i+2], a * b[i+3]
        sum2 *= _mm_load_ps(&b_row_vector[i+4]);  // a * b[i+4], a * b[i+5], a * b[i+6], a * b[i+7]
        result_row_vector[b_row_indices[i]]   += sum1[0];
        result_row_vector[b_row_indices[i+1]] += sum1[1];
        result_row_vector[b_row_indices[i+2]] += sum1[2];
        result_row_vector[b_row_indices[i+3]] += sum1[3];
        result_row_vector[b_row_indices[i+4]] += sum2[0];        
        result_row_vector[b_row_indices[i+5]] += sum2[1];        
        result_row_vector[b_row_indices[i+6]] += sum2[2];        
        result_row_vector[b_row_indices[i+7]] += sum2[3];        
    }
    for (; i < b_ellpack_cols; i++) // iterate through the last values (max 7 are left)
    {
        result_row_vector[b_row_indices[i]] += b_row_vector[i] * a_val;
    }
}


void sequential_multiplication(const void *a, const void *b, void *result) 
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
            scalar_multiplication_simd(a_val, b_row_vector, b_row_indices, result_matrix[a_row], b_matrix->ellpack_cols);
        }
    }
}
