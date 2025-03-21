#compiler flags
CFLAGS=-std=gnu11 -O3 -Wall -Wextra -Wpedantic

.PHONY: all clean debug msan

EXEC = matrix_multiplication

all: $(EXEC) 

$(EXEC): main.c V0/matr_mult_ellpack.c V1/matr_mult_ellpack_v1.c V2/matr_mult_ellpack_v2.c utils.c optimizations.c testing_functions.c
	$(CC) $(CFLAGS) -o $(EXEC) main.c V0/matr_mult_ellpack.c V1/matr_mult_ellpack_v1.c V2/matr_mult_ellpack_v2.c utils.c optimizations.c testing_functions.c

debug: main.c V0/matr_mult_ellpack.c V1/matr_mult_ellpack_v1.c V2/matr_mult_ellpack_v2.c utils.c optimizations.c testing_functions.c
	$(CC) $(CFLAGS) -g -o $(EXEC) main.c V0/matr_mult_ellpack.c V1/matr_mult_ellpack_v1.c V2/matr_mult_ellpack_v2.c utils.c optimizations.c testing_functions.c

asan: main.c V0/matr_mult_ellpack.c V1/matr_mult_ellpack_v1.c V2/matr_mult_ellpack_v2.c utils.c optimizations.c testing_functions.c
	$(CC) $(CFLAGS) -fsanitize=address -g -o $(EXEC) main.c V0/matr_mult_ellpack.c V1/matr_mult_ellpack_v1.c V2/matr_mult_ellpack_v2.c utils.c optimizations.c testing_functions.c
clean:
	@echo "Clean up"
	rm -f $(EXEC)