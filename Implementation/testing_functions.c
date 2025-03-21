#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include "V0/matr_mult_ellpack.h"
#include "V1/matr_mult_ellpack_v1.h"
#include "V2/matr_mult_ellpack_v2.h"
#include "utils.h"
#include <math.h>
#include <stdbool.h>
#include <unistd.h>

//testprogram for testing correctness of multiplication results of valid matrices


float random_float(float min, float max) {
    return min + ((float) rand() / (float) RAND_MAX) * (max - min);
}

//creates ellpack matrix with random float values between defined range
EllpackMatrix *create_random_ellpack_matrix(uint64_t rows, uint64_t cols, uint64_t ellpack_cols) {

    EllpackMatrix *matrix = allocate_ellpack_matrix(rows, cols, ellpack_cols);
    if(matrix == NULL){
        return NULL;
    }

    matrix->rows = rows;
    matrix->cols = cols;
    matrix->ellpack_cols = ellpack_cols;


    for (uint64_t i = 0; i < rows; i++) {
        int *used_indices = (int *) calloc(cols, sizeof(int));

            for (uint64_t j = 0; j < ellpack_cols; j++) {
            float random_f = random_float(-100,100);

            //after zero only zeros
            uint64_t k = j;
            if(random_f == 0.0f){
                while(k < ellpack_cols){
                    matrix->values[i][k] = 0.0f;
                    k++;
                }
                continue;
            }

            matrix->values[i][j] = random_float(-100,100);

            uint64_t rand_index;
            do {
                rand_index = rand() % cols;
            } while (used_indices[rand_index]); // unique col index

            matrix->indices[i][j] = rand_index;
            used_indices[rand_index] = 1;
        }
        free(used_indices);
    }

    return matrix;
}

float **allocate_2d_float_array(size_t rows, size_t cols) { //with calloc

    float **array = (float **) calloc(rows, sizeof(float *));
    if (array == NULL) {
        fprintf(stderr, ERR_MEMORY_ALLOCATION_FAILED_FAILED, "failed to allocate outer float array");
        return NULL;
    }

    for (size_t i = 0; i < rows; i++) {
        array[i] = (float *) calloc(cols, sizeof(float));
        if (array[i] == NULL) {
            fprintf(stderr, ERR_MEMORY_ALLOCATION_FAILED_FAILED, "failed to allocate inner arrays of float array");
            for (size_t j = 0; j < i; j++) {
                free(array[j]);
            }
            free(array);
            return NULL;
        }
    }

    return array;
}

void free_2d_float_array(float **array, size_t rows) {
    for (size_t i = 0; i < rows; ++i) {
        free(array[i]);
    }
    free(array);
}

//dump method for tests, writing multiple results to the same file
void test_dump_ellpack_matrix(FILE *file, EllpackMatrix *matrix) {

    // write <noRows>,<noCols>,<noEllpackCol> to the first line
    fprintf(file, "%"PRIu64",%"PRIu64",%"PRIu64"\n", matrix->rows, matrix->cols, matrix->ellpack_cols);

    // write matrix values to the second line
    for (uint64_t i = 0; i < matrix->rows; i++) {
        for (uint64_t j = 0; j < matrix->ellpack_cols; j++) {
            if (i == matrix->rows - 1 && j == matrix->ellpack_cols - 1) {
                fprintf(file, "%.1f", matrix->values[i][j]);
            } else {
                fprintf(file, "%.1f,", matrix->values[i][j]);
            }
        }
    }
    fprintf(file, "\n");

    // write indices to the third line
    for (uint64_t i = 0; i < matrix->rows; i++) {
        for (uint64_t j = 0; j < matrix->ellpack_cols; j++) {
            {
                if(i == matrix->rows - 1 && j == matrix->ellpack_cols - 1){
                    fprintf(file, "%"PRIu64"", matrix->indices[i][j]);
                }else{
                    fprintf(file, "%"PRIu64",", matrix->indices[i][j]);
                }
            }
        }
    }
    fprintf(file, "\n");

}

//dump method for tests, writing multiple results to the same file
char test_dump_result_to_ellpack(FILE* file, float **result_matrix, uint64_t rows, uint64_t cols)
{

    // check whether overflow occurred during multiplication
    for(uint64_t i = 0; i < rows; i++){
        for(uint64_t j = 0; j < cols; j++){
            if(isinf(result_matrix[i][j])){
                fprintf(stderr, ERR_OVERFLOW);
                return 'F';
            }
        }
    }
    // calculate ellpack_col
    uint64_t ellpack_cols = 0;
    for (uint64_t i = 0; i < rows; i++)
    {
        uint64_t nonZero = 0;
        for (uint64_t j = 0; j < cols; j++)
        {
            if (result_matrix[i][j] != 0)
            {
                nonZero++;
            }
        }
        if (nonZero > ellpack_cols)
        {
            ellpack_cols = nonZero;
        }
    }

    // write <noRows>,<noCols>,<noEllpackCol> to the first line
    fprintf(file, "%"PRIu64",%"PRIu64",%"PRIu64"\n", rows, cols, ellpack_cols);

    uint64_t* index_array = malloc(rows * ellpack_cols * sizeof(uint64_t));
    if (index_array == NULL) {
        fprintf(stderr, ERR_MEMORY_ALLOCATION_FAILED_FAILED, "failed to allocate memory for index_array");
        return 'F';
    }
    uint64_t array_index = 0;

       // write matrix values in ellpack-format to the second line
    for (uint64_t i = 0; i < rows; i++) {
        uint64_t pos = 0;
        for (uint64_t j = 0; j < cols; j++) {
            if (result_matrix[i][j] != 0) {
                if (pos > 0) {
                    fprintf(file, ",");
                }
                fprintf(file, "%.1f", result_matrix[i][j]);
                index_array[array_index++] = j;
                pos++;
            }
        }
        if (pos == ellpack_cols && i < rows - 1) //at the end of a row and not last row
        {
            fprintf(file, ",");
        }
        // fill unused places with *
        bool filled_unused = false;

        if (pos < ellpack_cols) {
            filled_unused = true;
        }

        while (pos < ellpack_cols) //row not full
        {
            if (pos > 0) {
                fprintf(file, ",*");
            } else {
                fprintf(file, "*");
            }
            index_array[array_index++] = -1;
            pos++;

        }
        if (filled_unused && i < rows - 1) { //at the end of a row but not last row
            fprintf(file, ",");
        }

    }
    fprintf(file, "\n");

    // write indices to the third line
    for (uint64_t i = 0; i < array_index; i++)
    {
        if (i > 0)
        {
            fprintf(file, ",");
        }
        if (index_array[i] == -1ULL) // comparison of integers of different signs: 'uint64_t'
        {
            fprintf(file, "*");
        }
        else
        {
            fprintf(file, "%"PRIu64"", index_array[i]);
        }
    }

    fprintf(file, "\n");
    free(index_array);
    return 'S';
}

//converts random gererated ellpack matrix to normal matrix for comparing with normal matrix multiplication
float** convert_ellpack_to_normal(EllpackMatrix* ellpackMatrix){
    uint64_t rows = ellpackMatrix->rows;
    uint64_t cols = ellpackMatrix->cols;
    uint64_t ellpack_cols = ellpackMatrix->ellpack_cols;
    float** value_array = ellpackMatrix->values;
    uint64_t ** indices_array = ellpackMatrix->indices;
    float** normal_matrix = allocate_2d_float_array(rows,cols);
    if(normal_matrix == NULL){
        exit(EXIT_FAILURE);
    }
    for(uint64_t i = 0; i < rows; i++){
        for(uint64_t j = 0; j < ellpack_cols; j++){
            normal_matrix[i][indices_array[i][j]] = value_array[i][j];
        }
    }
    return normal_matrix;
}

//standard matrix multiplication for comparison
float** normal_matrix_multiplication(float** matrix_a, float** matrix_b, uint64_t rows_a, uint64_t cols_a, uint64_t rows_b, uint64_t cols_b){
    float** result_matrix = allocate_2d_float_array(rows_a, cols_b);
    if(result_matrix == NULL || cols_a != rows_b){
        exit(EXIT_FAILURE);
    }
    for (uint64_t i= 0; i < rows_a; i++) {
        for (uint64_t j = 0; j < cols_b; j++) {
            for (uint64_t k = 0; k < cols_a; k++) {
                result_matrix[i][j] += matrix_a[i][k] * matrix_b[k][j];
            }
        }
    }
    return result_matrix;
}

//comparing results between different multiplications
bool compare(float** ellpack_result, float** normal_result,uint64_t rows, uint64_t cols){
    for (uint64_t i = 0; i < rows; i++) {
        for (uint64_t j = 0; j < cols; j++) {
            if (fabs(ellpack_result[i][j] - normal_result[i][j]) > 1.0) {
                fprintf(stderr, "should: %f, but is: %f\n", normal_result[i][j], ellpack_result[i][j]);
                return false;
            }
        }
    }
    return true;
}

//testing multiplication results with normal matrix multiplication results
double run_multiplication_test(FILE * file, uint64_t rows_a, uint64_t cols_a, uint64_t ellpack_cols_a, uint64_t rows_b, uint64_t cols_b, uint64_t ellpack_cols_b, int version){
    fprintf(file, "testing matrix multiplication with version %i:\n", version);
    fprintf(file, "generating random matrices of float values between -100 and 100: \n");
    fprintf(file, "\n");

    struct timespec start, end;
    double elapsed_time = 0;
    double elapsed_time_normal = 0;

    EllpackMatrix* test_a = create_random_ellpack_matrix(rows_a,cols_a,ellpack_cols_a);
    EllpackMatrix* test_b = create_random_ellpack_matrix(rows_b,cols_b,ellpack_cols_b);
    if(test_a == NULL || test_b == NULL){
        exit(EXIT_FAILURE);
    }
    float** result_pointer = allocate_2d_float_array(rows_a,cols_b);
    if(result_pointer == NULL){
        exit(EXIT_FAILURE);
    }
    if (version == 0) 
    {
        clock_gettime(CLOCK_MONOTONIC, &start);
        matr_mult_ellpack(test_a, test_b, result_pointer);
        sleep(1);
        clock_gettime(CLOCK_MONOTONIC, &end);
        elapsed_time = (end.tv_sec - start.tv_sec) * 1.0e9 + (end.tv_nsec - start.tv_nsec);
    }
    else if (version == 1) 
    {
        clock_gettime(CLOCK_MONOTONIC, &start);
        matr_mult_ellpack_v1(test_a, test_b, result_pointer);
        sleep(1);
        clock_gettime(CLOCK_MONOTONIC, &end);
        elapsed_time = (end.tv_sec - start.tv_sec) * 1.0e9 + (end.tv_nsec - start.tv_nsec);
    }
    else 
    {
        clock_gettime(CLOCK_MONOTONIC, &start);
        matr_mult_ellpack_v2(test_a, test_b, result_pointer);
        sleep(1);
        clock_gettime(CLOCK_MONOTONIC, &end);
        elapsed_time = (end.tv_sec - start.tv_sec) * 1.0e9 + (end.tv_nsec - start.tv_nsec);
    }

    fprintf(file, "matrix_a: \n");
    test_dump_ellpack_matrix(file, test_a);
    fprintf(file, "\n");
    fprintf(file, "matrix_b: \n");
    test_dump_ellpack_matrix(file, test_b);
    fprintf(file, "\n");
    fprintf(file, "result: \n");
    test_dump_result_to_ellpack(file, result_pointer, rows_a,cols_b);
    fprintf(file, "\n");

    //test multiplication result
    float** normal_a = convert_ellpack_to_normal(test_a);
    float** normal_b = convert_ellpack_to_normal(test_b);

    //normal multiplication performance
    clock_gettime(CLOCK_MONOTONIC, &start);
    float** normal_res = normal_matrix_multiplication(normal_a,normal_b,rows_a,cols_a,rows_b,cols_b);
    sleep(1);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time_normal = (end.tv_sec - start.tv_sec) * 1.0e9 + (end.tv_nsec - start.tv_nsec);

    bool test_res = compare(result_pointer, normal_res, rows_a,cols_b);
    if(test_res){
        fprintf(file, "multiplication successful\n");
        fprintf(file, "accepted float variation tolerance was: 1.0\n");
    }else{
        fprintf(file, "multiplication failed\n");
        fprintf(file, "accepted float variation tolerance was: 1.0, variations are shown on the console, adjust the tolerance\n");
    }

    fprintf(file, "benchmarking: \n");
    fprintf(file, "ellpack matrix multiplication took: %f seconds\n", elapsed_time / 1.0e9 - 1);
    fprintf(file, "normal matrix multiplication took: %f seconds\n", elapsed_time_normal / 1.0e9 - 1);



    free_2d_float_array(result_pointer,rows_a);
    free_2d_float_array(normal_a,rows_a);
    free_2d_float_array(normal_b,rows_b);
    free_2d_float_array(normal_res,rows_a);
    free_ellpack_matrix(test_a);
    free_ellpack_matrix(test_b);

    return elapsed_time / 1.0e9 - 1;
}

//method invoked by main.c
int execute_tests(int version){

    srand(time(NULL));
    double total_time = 0;

    char* filename;
    if (version == 0) 
    {
        filename = "test_v0.txt";
    }
    else if (version == 1)
    {
        filename = "test_v1.txt";
    }
    else 
    {
        filename = "test_v2.txt";
    }
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        fprintf(stderr, ERR_OPEN_FILE_FAILED, filename);
        return -1;
    }
    int i = 0;

    //multiplication test + performance
    for(; i < 3; i++){
        fprintf(file, "round: %i\n", i);
        total_time += run_multiplication_test(file,20,20,10,20,20,10,version);
        fprintf(file, "\n\n\n");
    }

    double average_execution_time = total_time / i;
    fprintf(file, "Average time of ellpack-mul: %f seconds\n", average_execution_time);

    fclose(file);
    return 0;
}