#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "utils.h"
#include "V0/matr_mult_ellpack.h"
#include "V1/matr_mult_ellpack_v1.h"
#include "V2/matr_mult_ellpack_v2.h"
#include "testing_functions.h"


static struct option long_options[] = {
    {"iterations", optional_argument, 0, 'B'},
    {"version", required_argument, 0, 'V'},
    {"matrix_a", required_argument, 0, 'a'},
    {"matrix_b", required_argument, 0, 'b'},
    {"output", required_argument, 0, 'o'},
    {"help", no_argument, 0, 'h'},
    {"test", no_argument, 0, 't'},
    {0, 0, 0, 0}};

void print_usage(void)
{
    printf("Here is the command usage documentation: ");
    printf("./matrix_multiplication --matrix_a <file> --matrix_b <file> --output <file> --version [0-2]\n");
    printf("options:\n");
    printf("  -B, --iterations=[N]     The number of iterations for benchmarking\n");
    printf("  -V, --version [0-2]  The implementation version\n");
    printf("  -a, --matrix_a <file>    Input file for matrix A\n");
    printf("  -b, --matrix_b <file>    Input file for matrix B\n");
    printf("  -o, --output <file>      Output file for the result matrix\n");
    printf("  -t, --test               Run multiplication tests\n");
    printf("  -h, --help               Display this help message\n");
}

int main(int argc, char **argv)
{
    int opt;
    int option_index = 0;
    char *a_filename = NULL;
    char *b_filename = NULL;
    char *output_filename = NULL;
    unsigned int version = 0;    // use implementation 0 by default (naive algorithm)
    unsigned int iterations = 1; // iterate only once by default

    if (strcmp(argv[0], "./matrix_multiplication") != 0)
    {
        fprintf(stderr, "Error: The program name entered is incorrect.\n");
        print_usage();
        exit(EXIT_FAILURE);
    }

    // parse the options
    while ((opt = getopt_long(argc, argv, "B::V:a:b:o:h:t", long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 'V':
        {
            char *v_endptr;
            version = strtol(optarg, &v_endptr, 10);
            break;
        }
        case 'B':
        {
            char *i_endptr;
            iterations = strtol(optarg, &i_endptr, 10);
            break;
        }
        case 'a':
            a_filename = optarg;
            break;
        case 'b':
            b_filename = optarg;
            break;
        case 'o':
            output_filename = optarg;
            break;
        case 'h':
            print_usage();
            exit(EXIT_SUCCESS);
        case 't':
            execute_tests(0);
            execute_tests(1);
            execute_tests(2);
            exit(EXIT_SUCCESS);
        default:
            fprintf(stderr, "Error: Invalid argument detected.\n");
            print_usage();
            exit(EXIT_FAILURE);
        }
    }

    if (!a_filename || !b_filename || !output_filename)
    {
        fprintf(stderr, "Error: Missing required arguments.\n");
        exit(EXIT_FAILURE);
    }

    // run the matrix multiplication for a certain amount of iterations
    EllpackMatrix *m1 = load_ellpack_matrix(a_filename);
    if (m1 == NULL) 
    {
        fprintf(stderr, "Error: Failed to allocate m1 matrix\n");
        exit(EXIT_FAILURE);
    }
    EllpackMatrix *m2 = load_ellpack_matrix(b_filename);
    if (m2 == NULL)
    {
        free_ellpack_matrix(m1);
        fprintf(stderr, "Error: Failed to allocate m2 matrix\n");
        exit(EXIT_FAILURE);
    }


    struct timespec start, end;

    unsigned int counter = 0;
    double elapsed_time = 0;
    while (counter < iterations)
    {
        if (version > 2) 
        {
            fprintf(stderr, "Error: The version number \"%d\" is invalid.\n", version);
            free_ellpack_matrix(m1);
            free_ellpack_matrix(m2);
            exit(EXIT_FAILURE);
        }

        float **res = allocate_matrix_array(m1->rows, m2->cols);
        if (res == NULL)
        {
            fprintf(stderr, "Error: Failed to allocate res.\n");
            free_ellpack_matrix(m1);
            free_ellpack_matrix(m2);
            exit(EXIT_FAILURE);
        }
        if (version == 0)
        {
            clock_gettime(CLOCK_MONOTONIC, &start);
            matr_mult_ellpack(m1, m2, res);
            sleep(1);
            clock_gettime(CLOCK_MONOTONIC, &end);
        }
        else if (version == 1)
        {
            clock_gettime(CLOCK_MONOTONIC, &start);
            matr_mult_ellpack_v1(m1, m2, res);
            sleep(1);
            clock_gettime(CLOCK_MONOTONIC, &end);
        }
        else if (version == 2)
        {
            clock_gettime(CLOCK_MONOTONIC, &start);
            matr_mult_ellpack_v2(m1, m2, res);
            sleep(1);
            clock_gettime(CLOCK_MONOTONIC, &end);
        }
        char dumped = dump_result_to_ellpack(output_filename, res, m1->rows, m2->cols);
        free_matrix_array(m1->rows, res);
        if (dumped == 'F')
        {
            free_ellpack_matrix(m1);
            free_ellpack_matrix(m2);
            exit(EXIT_FAILURE);
        }
        counter++;
        elapsed_time += (end.tv_sec - start.tv_sec) * 1.0e9 + (end.tv_nsec - start.tv_nsec);
    }

    // calculate the duration of matr_mult_ellpack after n iterations
    printf("Version %d average elapsed time per iteration: %f seconds (%d iterations ran)\n", version, (elapsed_time / iterations / 1.0e9) - 1.0, iterations);

    free_ellpack_matrix(m1);
    free_ellpack_matrix(m2);

    exit(EXIT_SUCCESS);
}
