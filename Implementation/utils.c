#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <limits.h>
#include "utils.h"
#include <math.h>
#include <stdbool.h>

EllpackMatrix *load_ellpack_matrix(const char *filename)
{
    FILE *file;
    char *line = NULL;
    size_t len = 0;

    if (NULL == (file = fopen(filename, "r")))
    {
        fprintf(stderr, ERR_OPEN_FILE_FAILED, filename);
        return NULL;
    }

    // ######### FIRST SECTION #########

    // process first line
    // format: "<rows>,<cols>,<noEllpackRows>"
    if (getline(&line, &len, file) == -1)
    {
        fprintf(stderr, ERR_READ_LINE_FAILED, 1, filename);
        free(line);
        fclose(file);
        return NULL;
    }

    // split the line into tokens by ","
    // only process the first 3 tokens and convert them into uint64
    char *token = strtok(line, ",");
    uint64_t i = 0;
    uint64_t matr_params[3];
    matr_params[0] = -1; // dummy val
    matr_params[1] = -1;
    matr_params[2] = -1; 
    while (token != NULL && i < 3)
    {
        char cleaned_token[64] = {0};
        remove_invalid_chars(cleaned_token, token);

        if (convert_to_uint64(cleaned_token, &matr_params[i]) != 'S')
        {
            fprintf(stderr, ERR_CONVERT_UINT64_FAILED, token);
            free(line);
            fclose(file);
            return NULL;
        }
        // the dimensions of the matrix cannot be zero
        if (matr_params[i] == 0)
        {
            fprintf(stderr, ERR_ZERO_DIMENSION);
            free(line);
            fclose(file);
            return NULL;
        }

        i++;
        token = strtok(NULL, ",");
    }

    if (matr_params[2] > matr_params[1]) {
        fprintf(stderr, ERR_INVALID_ELLPACK_COLS);
        free(line);
        fclose(file);
        return NULL;
    }

    // missing tokens in the first line
    if (i != 3)
    {
        fprintf(stderr, ERR_UNEXPECTED_TOKEN_NUMBER, (uint64_t)3, i);
        free(line);
        fclose(file);
        return NULL;
    }

    // initialize the ellpack matrix
    EllpackMatrix *matrix = allocate_ellpack_matrix(matr_params[0], matr_params[1], matr_params[2]);
    if (matrix == NULL)
    {
        free(line);
        fclose(file);
        return NULL;
    }

    // ######### SECOND SECTION #########

    // read the second line, which contains the values of the ELLPACK matrix
    if (getline(&line, &len, file) == -1)
    {
        fprintf(stderr, ERR_READ_LINE_FAILED, 2, filename);
        free_ellpack_matrix(matrix);
        free(line);
        fclose(file);
        return NULL;
    }

    // split the line into tokens by ","
    token = strtok(line, ",");
    uint64_t row = 0;
    uint64_t col = 0;
    while (token != NULL)
    {
        if (row == matrix->rows)
        {
            // stop reading if all rows are read already
            break;
        }

        char cleaned_token[64] = {0};
        remove_invalid_chars(cleaned_token, token);
        if (strcmp(cleaned_token, "*") == 0)
        {
            matrix->values[row][col] = 0.0F; // default '*' to 0.0
        }
        else
        {
            char *endptr;
            float value = strtof(cleaned_token, &endptr);
            if (*endptr != '\0')
            {
                fprintf(stderr, ERR_CONVERT_TO_FLOAT_FAILED, token);
                free_ellpack_matrix(matrix);
                free(line);
                fclose(file);
                return NULL;
            }
            if (isinf(value)) {
                fprintf(stderr, ERR_CONVERT_TO_FLOAT_FAILED, token);
                free_ellpack_matrix(matrix);
                free(line);
                fclose(file);
                return NULL;
            }
            if (isnan(value)) {
                fprintf(stderr, ERR_CONVERT_TO_FLOAT_FAILED, token);
                free_ellpack_matrix(matrix);
                free(line);
                fclose(file);
                return NULL;
            }
            matrix->values[row][col] = value;
        }
        col++;

        // reset the col index and update row index with respect to `ellpack_cols`
        if (col == matrix->ellpack_cols)
        {
            col = 0;
            row++;
        }
        token = strtok(NULL, ",");
    }

    // check if the number of tokens is correct
    if (row != matrix->rows)
    {
        fprintf(stderr, ERR_UNEXPECTED_ROW_NUMBER, matrix->rows, row);
        free_ellpack_matrix(matrix);
        free(line);
        fclose(file);
        return NULL;
    }

    // ######### THIRD SECTION #########

    // read the third line, which contains the indices of the ELLPACK matrix
    if (getline(&line, &len, file) == -1)
    {
        fprintf(stderr, ERR_READ_LINE_FAILED, 3, filename);
        free_ellpack_matrix(matrix);
        free(line);
        fclose(file);
        return NULL;
    }

    // split the line into tokens by ","
    token = strtok(line, ",");
    row = 0;
    col = 0;
    char *appeared = (char *)calloc(matrix->cols, sizeof(char));
    if (appeared == NULL) {
        fprintf(stderr, "Failed to allocate `appeared` array for checking duplicate indices.");
        free_ellpack_matrix(matrix);
        free(line);
        fclose(file);
        return NULL;
    }
    while (token != NULL)
    {
        if (row == matrix->rows)
        {
            // stop reading if all rows are read already
            break;
        }

        if (strcmp(token, "*") == 0)
        {
            matrix->indices[row][col] = 0; // default '*' to 0
        }
        else
        {
            char cleaned_token[64] = {0};
            remove_invalid_chars(cleaned_token, token);
            uint64_t value = 0;
            if (convert_to_uint64(cleaned_token, &value) != 'S')
            {
                fprintf(stderr, ERR_CONVERT_UINT64_FAILED, token);
                free_ellpack_matrix(matrix);
                free(line);
                fclose(file);
                free(appeared);
                return NULL;
            }
            if (value >= matrix->cols) {
                fprintf(stderr, ERR_INVALID_INDEX, value);
                free_ellpack_matrix(matrix);
                free(line);
                fclose(file);
                free(appeared);
                return NULL;
            }
            // check for duplicate indices
            if (appeared[value] == 1 && matrix->values[row][col] != 0.0F) {
                fprintf(stderr, ERR_DUPLICATE_INDEX, value);
                free(appeared);
                free_ellpack_matrix(matrix);
                free(line);
                fclose(file);
                return NULL;
            }
            matrix->indices[row][col] = value;
            appeared[value] = 1;
        }
        col++;
        if (col == matrix->ellpack_cols)
        {
            col = 0;
            row++;
            free(appeared);
            appeared = (char *)calloc(matrix->cols, sizeof(char));
            if (appeared == NULL) {
                fprintf(stderr, "Failed to allocate `appeared` array for checking duplicate indices.");
                free_ellpack_matrix(matrix);
                free(line);
                fclose(file);
                return NULL;
            }
        }
        token = strtok(NULL, ",");
    }

    // free the array for checking duplicate indices
    free(appeared);

    // check if the number of tokens is correct
    if (row != matrix->rows)
    {
        fprintf(stderr, ERR_UNEXPECTED_ROW_NUMBER, matrix->rows, row);
        free_ellpack_matrix(matrix);
        free(line);
        fclose(file);
        return NULL;
    }

    free(line);
    fclose(file);
    return matrix;
}

//convert normal result matrix to ellpack format and write it to result.txt
char dump_result_to_ellpack(const char *filename, float **result_matrix, uint64_t rows, uint64_t cols) {

    // check whether overflow occurred during multiplication
    for (uint64_t i = 0; i < rows; i++) {
        for (uint64_t j = 0; j < cols; j++) {
            if (isinf(result_matrix[i][j])) {
                fprintf(stderr, ERR_OVERFLOW);
                return 'F';
            }
        }
    }

    // calculate ellpack_cols
    uint64_t ellpack_cols = 0;
    for (uint64_t i = 0; i < rows; i++) {
        uint64_t nonZero = 0;
        for (uint64_t j = 0; j < cols; j++) {
            if (result_matrix[i][j] != 0) {
                nonZero++;
            }
        }
        if (nonZero > ellpack_cols) {
            ellpack_cols = nonZero;
        }
    }

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, ERR_OPEN_FILE_FAILED, filename);
        return 'F';
    }

    // write <noRows>,<noCols>,<noEllpackCol> to the first line
    fprintf(file, "%"PRIu64",%"PRIu64",%"PRIu64"\n", rows, cols, ellpack_cols);

    uint64_t *index_array = malloc(rows * ellpack_cols * sizeof(uint64_t));
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
    fclose(file);

    free(index_array);
    return 'S';
}

//write ellpack matrix directly to output file
void dump_ellpack_matrix(const char *filename, const EllpackMatrix *matrix)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        fprintf(stderr, ERR_OPEN_FILE_FAILED, filename);
        return;
    }

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

    fclose(file);
}

/*
 * Remove white spaces and `\n` from a string
 */
void remove_invalid_chars(char *dest, const char *src)
{
    while (*src)
    {
        if (strcmp((const char *)src, " ") != 0 && strcmp((const char *)src, "\n") != 0)
        {
            *dest++ = *src;
        }
        src++;
    }
    *dest = '\0';
}

/*
 * Helper method to convert a string into uint64
 */
char convert_to_uint64(const char *str, uint64_t *value)
{
    // check for negative sign
    for (u_int64_t i = 0; i < sizeof(str); i++) {
        if (str[i] == '-') {
            return 'F';
        }
    }

    // the string is empty
    if (str[0] == '\0')
    {
        return 'F';
    }

    char *endptr;
    errno = 0;
    uint64_t ull_value = strtoull(str, &endptr, 10);

    // input is invalid
    if (*endptr != '\0')
    {
        return 'F';
    }
    // overflow during the conversion
    else if (errno == ERANGE && ull_value == ULLONG_MAX)
    {
        return 'F';
    }

    *value = (uint64_t)ull_value;
    return 'S';
}

/*
 * Helper method to allocate memory for an ELLPACK matrix
 * Using calloc() to initialize the values and indices array so the values default to zero
 */
EllpackMatrix *allocate_ellpack_matrix(uint64_t rows, uint64_t cols, uint64_t ellpack_cols)
{
    // allocate memory for the EllpackMatrix structure
    EllpackMatrix *matrix = (EllpackMatrix *)malloc(sizeof(EllpackMatrix));
    if (matrix == NULL)
    {
        fprintf(stderr, ERR_MEMORY_ALLOCATION_FAILED_FAILED, "matrix structure");
        return NULL;
    }

    // initialize the dimensions
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->ellpack_cols = ellpack_cols;

    // allocate memory for the values and indices arrays
    matrix->values = (float **)calloc(rows, sizeof(float *));
    matrix->indices = (uint64_t **)calloc(rows, sizeof(uint64_t *));
    if (matrix->values == NULL || matrix->indices == NULL)
    {
        fprintf(stderr, ERR_MEMORY_ALLOCATION_FAILED_FAILED, "matrix values/indices array");
        free(matrix->values);
        free(matrix->indices);
        free(matrix);
        return NULL;
    }

    // allocate memory for each row in the values and indices arrays
    for (uint64_t i = 0; i < rows; i++)
    {
        matrix->values[i] = (float *)calloc(ellpack_cols, sizeof(float));
        matrix->indices[i] = (uint64_t *)calloc(ellpack_cols, sizeof(uint64_t));
        if (matrix->values[i] == NULL || matrix->indices[i] == NULL)
        {
            fprintf(stderr, ERR_MEMORY_ALLOCATION_FAILED_FAILED, "inner array of matrix values/indices array");
            for (uint64_t j = 0; j <= i; j++)
            {
                free(matrix->values[j]);
                free(matrix->indices[j]);
            }
            free(matrix->values);
            free(matrix->indices);
            free(matrix);
            return NULL;
        }
    }

    return matrix;
}

/*
 * Helper method to free memory of an ELLPACK matrix
 */
void free_ellpack_matrix(EllpackMatrix *matrix)
{
    if (matrix != NULL)
    {
        for (uint64_t i = 0; i < matrix->rows; i++)
        {
            free(matrix->values[i]);
            free(matrix->indices[i]);
        }
        free(matrix->values);
        free(matrix->indices);
        free(matrix);
    }
}

float **allocate_matrix_array(uint64_t rows, uint64_t cols)
{
    float **arr = (float **)calloc(rows, sizeof(float *));
    if (arr == NULL)
    {
        fprintf(stderr, "Memory allocation failed for row pointers\n");
        return NULL;
    }

    for (uint64_t i = 0; i < rows; i++)
    {
        arr[i] = (float *)calloc(cols, sizeof(float));
        if (arr[i] == NULL)
        {
            fprintf(stderr, "Memory allocation failed for columns in row %"PRIu64"\n", i);
            // free previously allocated memory before exiting
            for (uint64_t j = 0; j < i; j++)
            {
                free(arr[j]);
            }
            free(arr);
            return NULL;
        }
    }

    return arr;
}

void free_matrix_array(uint64_t rows, float **arr)
{
    for (uint64_t i = 0; i < rows; i++)
    {
        free(arr[i]);
        arr[i] = NULL;
    }
    free(arr);
    arr = NULL;
}
