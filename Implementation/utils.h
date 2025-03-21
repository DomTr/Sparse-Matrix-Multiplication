#ifndef FINAL_UTILS_H
#define FINAL_UTILS_H
#include <inttypes.h>
// Error message format strings

#define ERR_OPEN_FILE_FAILED "Error: Failed to open %s\n"
#define ERR_READ_LINE_FAILED "Error: Failed to read the line no.%d in %s\n"
#define ERR_CONVERT_UINT64_FAILED "Error: Failed to convert '%s' to uint64_t\n"
#define ERR_CONVERT_TO_FLOAT_FAILED "Error: Could not convert token '%s' to float\n"
#define ERR_MEMORY_ALLOCATION_FAILED_FAILED "Error: Memory allocation failed for %s, the memory needed is possibly too large\n"
#define ERR_UNEXPECTED_TOKEN_NUMBER "Error: Expected %"PRIu64" tokens but found %"PRIu64"\n"
#define ERR_INVALID_MATRIX_DIMENSIONS "Error: columns in matrix A (%" PRIu64 ") != rows in matrix B (%"PRIu64")\n"
#define ERR_ZERO_DIMENSION "Error: the dimensions of the matrix cannot be 0"
#define ERR_UNEXPECTED_ROW_NUMBER "Error: Expected %"PRIu64" rows but %"PRIu64" rows were read"
#define ERR_OVERFLOW "Error: Overflow occurred during multiplication\n"
#define ERR_INVALID_INDEX "Error: The index %"PRIu64" is out of bounds\n"
#define ERR_INVALID_ELLPACK_COLS "Error: EllpackCols should be less than the column dimension of the matrix\n"
#define ERR_DUPLICATE_INDEX "Error: The index %"PRIu64" is duplicated\n"

// EllpackMatrix structure

typedef struct
{
    uint64_t rows;
    uint64_t cols;
    uint64_t ellpack_cols;
    float **values;
    uint64_t **indices;
} EllpackMatrix;

// Helper methods

void remove_invalid_chars(char *dest, const char *src);

char convert_to_uint64(const char *str, uint64_t *value);

EllpackMatrix *load_ellpack_matrix(const char *filename);

char dump_result_to_ellpack(const char *filename, float **result_matrix, uint64_t rows, uint64_t cols);

void dump_ellpack_matrix(const char *filename, const EllpackMatrix *matrix);

EllpackMatrix *allocate_ellpack_matrix(uint64_t rows, uint64_t cols, uint64_t ellpack_cols);

void free_ellpack_matrix(EllpackMatrix *matrix);

float **allocate_matrix_array(uint64_t rows, uint64_t cols);

void free_matrix_array(uint64_t rows, float **arr);

#endif
