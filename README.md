# Sparse_Ellpack_Matrix_Mult
Sparse Matrix Multiplication in ELLPACK format, a university project of 3 people

How to run the program:
1. Navigate to Implementation directory. If you are on intel processor, skip steps 2 and 3, just run 'make' inside terminal. matrix_multiplication executable will be created. Docker is required for ARM-based processors because SIMD instructions are used.
2. sudo docker build --platform linux/amd64 -t make-image .
This will call make on the whole program
Then: 
3. docker run --platform linux/amd64 -it make-image bash 
4. Run ` ./matrix_multiplication --matrix_a InputData/matrix_100_100_25_0_a --matrix_b InputData/matrix_100_100_25_0_b --output result.txt -V2`
* matrix_multiplication is name of the executable 
* --matrix_a is command line argument for first matrix
* --matrix_b is command line for second matrix. Arguments must be files where matrices are stored in ELLPACK format.
* --output 'file_name' specifies where to store the product of the multiplied matrices. The resulting matrix is in also in ELLPACK.
* -V2 means running with optimization level 2. There are 3 levels, -V1 - no optimization, -V2 - optimization with SIMD, -V3 - optimization with SIMD and threads.
5. For users running on Docker, you can type ls in bash, then you will see your results file (here it is results.txt). Type cat <res_file_name> and you will see the result.
