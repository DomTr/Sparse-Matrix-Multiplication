#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "matr_mult_ellpack_v2.h"

#define NUM_THREADS 5

/*
 * Main function to multiply two EllpackMatrix and save it to the result pointer
 */
void matr_mult_ellpack_v2(const void *a, const void *b, void *result)
{
    EllpackMatrix *a_matrix = (EllpackMatrix *)a;
    if (a_matrix->rows <= 5 * NUM_THREADS)
    {
        sequential_multiplication(a, b, result);
    }
    else
    {
        parallel_multiplication(a, b, result);
    }
}

/*
 * Function to multiply matrices using threads.
 */
void parallel_multiplication(const void *a, const void *b, void *result)
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

    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    uint64_t chunk_size = a_matrix->rows / NUM_THREADS;

    // thread creation
    for (uint64_t i = 0; i < NUM_THREADS; i++)
    {
        thread_data[i].start = i * chunk_size;
        thread_data[i].end = (i == NUM_THREADS - 1) ? a_matrix->rows : (i + 1) * chunk_size; // the last thread can maximum go to the end
        thread_data[i].a_matrix = a_matrix;
        thread_data[i].b_matrix = b_matrix;
        thread_data[i].result = result_matrix;
        pthread_create(&threads[i], NULL, matr_mult_one_thread, (void *)&thread_data[i]);
    }

    for (int i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], NULL);
    }
}

/*
 * matr_mult_one_thread is a function that gets executed by a single thread. arg: ThreadData.
 * main difference from sequential multiplication is that a_row starts starts not from 0 and ends at a_rows, but from data->start and ends before data->end.
 */
void *matr_mult_one_thread(void *arg)
{
    ThreadData *data = (ThreadData *)arg;
    EllpackMatrix *a_matrix = data->a_matrix;
    EllpackMatrix *b_matrix = data->b_matrix;
    float **result_matrix = data->result;
    for (uint64_t a_row = data->start; a_row < data->end; a_row++)
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
    return NULL;
}
