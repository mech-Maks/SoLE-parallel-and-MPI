#include "matrix_operations/matrix_operations.h"

int main(int argc, char const *argv[]) {
/* args: matrix_size, block_size, threads_amount, print_size, formula_type, (optional) filename */
    feenableexcept (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW );
    
    double *A_matrix = NULL;

    double *B_vector = NULL;
    double *X_vector = NULL;
    double *X_true = NULL;

    int *X_indeces = NULL;

    int matrix_size = 0;
    int block_size = 0;
    int print_size = 0;
    int formula_type = 0;
    int i = 0;
    int threads_amount = 0;

    double time_res = 0;
    double time_discrepency = 0;
    double temp_double = 0;
    double residual = 0;

    pthread_t *tids;
    arguments *args;
    mainElem *main_elems;
    pthread_barrier_t barrier;

    clock_t start, end;

    if (argc < 6) {
        printf("Not enough argumenst\n");
        return -1;
    }

    if (argc > 7) {
        printf("Too many arguments\n");
        return -1;
    }

    if (sscanf(argv[1], "%d", &matrix_size) != 1    ||
        sscanf(argv[2], "%d", &block_size) != 1     ||
        sscanf(argv[3], "%d", &threads_amount) != 1 ||
        sscanf(argv[4], "%d", &print_size) != 1     ||
        sscanf(argv[5], "%d", &formula_type) != 1)
    {
        printf("Put numbers, not strings\n");
        return -1;
    }

    if (formula_type == 0 && argc < 6) {
        printf("No input-file\n");
        return -1;
    }

    if (matrix_size <= 0 ||
        block_size <= 0  ||
        print_size < 0   ||
        threads_amount <= 0  ) 
    {
        printf("Put in possitive numbers\n");
        return -1;
    }

    if (matrix_size < block_size) block_size = matrix_size;
    if (matrix_size < print_size) print_size = matrix_size;
    if (matrix_size < threads_amount) threads_amount = matrix_size;

    if (formula_type < 0 || formula_type > 4) {
        printf("Unknown formula type\n");
        return -1;
    }

    A_matrix = (double*)malloc(matrix_size*matrix_size*sizeof(double));
    if (!A_matrix) {
        printf("Компьютер умер\n");
        return -1;
    }

    pthread_barrier_init(&barrier, 0, threads_amount);

    X_indeces = (int*)calloc(matrix_size, sizeof(int));
    X_vector = (double*)calloc(matrix_size, sizeof(double));
    X_true = (double*)calloc(matrix_size, sizeof(double));

    B_vector = (double*)malloc(matrix_size*sizeof(double));

    for (i=0; i<matrix_size; i++) {
        X_indeces[i] = i;

        if (i%2 == 0) {
            X_true[i] = 1;
        } else {
            X_true[i] = 0;
        }
    }
    
    args = new arguments[threads_amount];
    tids = (pthread_t*)malloc(threads_amount*sizeof(pthread_t));
    main_elems = new mainElem[threads_amount];

    for (i=0; i<threads_amount; i++) {
        args[i].A = A_matrix;
        args[i].B = B_vector;
        args[i].X_indeces = X_indeces;
        args[i].formula_type = formula_type;
        args[i].matrix_size = matrix_size;
        args[i].block_size = block_size;
        args[i].t_index = i;
        args[i].t_amount = threads_amount;
        args[i].barrier = &barrier;
        args[i].main_elems = main_elems;
        args[i].method_does_not_work = 1;

        main_elems[i].value = 0.0;
        main_elems[i].index = 0;
        main_elems[i].error = 1;
    }

    if (formula_type == 0) {
        if ( !buildMatrixFromFile(A_matrix, matrix_size, formula_type, argv[6], B_vector) ) {
            printf("Failed while building matrix\n");
            if (A_matrix) free(A_matrix);
            if (B_vector) free(B_vector);
            if (X_vector) free(X_vector);
            if (X_indeces) free(X_indeces);
            if (X_true) free(X_true);
            pthread_barrier_destroy(&barrier);
            delete[] args;
            delete[] main_elems;
            free(tids);
            
            return -1;
        }
    }

    for (i=0; i<threads_amount-1; i++) {
        if (pthread_create(tids+i, 0, &thread_func, args+i)) {
            printf("Cannot create thread\n");
            if (A_matrix) free(A_matrix);
            if (B_vector) free(B_vector);
            if (X_vector) free(X_vector);
            if (X_indeces) free(X_indeces);
            if (X_true) free(X_true);
            pthread_barrier_destroy(&barrier);
            delete[] args;
            delete[] main_elems;
            free(tids);

            return -1;
        }
    }
    thread_func(args+threads_amount-1);

    for (int t=0; t<threads_amount-1; t++) {
        pthread_join(tids[t], 0);
    }

    for (i=0; i<threads_amount; i++) {
        if (args[i].method_does_not_work < 0) {
            printf("Найдена необратимая строка, метод гаусса неприменим\n");
            if (A_matrix) free(A_matrix);
            if (B_vector) free(B_vector);
            if (X_vector) free(X_vector);
            if (X_indeces) free(X_indeces);
            if (X_true) free(X_true);
            pthread_barrier_destroy(&barrier);
            delete[] args;
            delete[] main_elems;
            free(tids);


            return -1;
        }
    }
    
    time_res = args[0].elapsed_total;

    // printMatrix(A_matrix, matrix_size, matrix_size, print_size);
    

    for (i=0; i<matrix_size; i++) {
        X_vector[ X_indeces[i] ] = B_vector[i];
    }

    // printMatrix(X_vector, matrix_size, 1, print_size);


    // Восстанавливаем матрицу A
    if (formula_type == 0) {
        if ( !buildMatrixFromFile(A_matrix, matrix_size, formula_type, argv[6], B_vector) ) {
            printf("Failed while building matrix\n");
            if (A_matrix) free(A_matrix);
            if (B_vector) free(B_vector);
            if (X_vector) free(X_vector);
            if (X_indeces) free(X_indeces);
            if (X_true) free(X_true);
            pthread_barrier_destroy(&barrier);
            delete[] args;
            delete[] main_elems;
            free(tids);


            return -1;
        }
    } 
    else {
        for (i=0; i<matrix_size; i++) {
            for (int j=0; j<matrix_size; j++) {
                A_matrix[i*matrix_size + j] = buildMatrix(formula_type, matrix_size, i, j);
            }

            B_vector[i] = sum(A_matrix, matrix_size, i);
        }
    }

    if (threads_amount == 1 && matrix_size > 50) {
        residual = -1;
        time_discrepency = 0;
    }
    else {
        start = clock();
        residual = discrepancy(A_matrix, X_vector, B_vector, matrix_size);
        end = clock();

        time_discrepency = end-start;
        time_discrepency = (double)time_discrepency/CLOCKS_PER_SEC;
    }

    subtractMatrices(X_true, X_vector, matrix_size, 1);

    temp_double = matrixNorm(X_true, matrix_size, 1);
    printf("final norm: %10.3e\n", temp_double);

    printf("%.2f ms потребовалось на поиск невязки\n\n", time_discrepency);

    // printf("\n\n==========================================\n\n");
    // for (i=0; i<threads_amount; i++) {
    //     printf("%d\telapsed_thread: %lf\n", i, args[i].elapsed_thread);
    //     printf("\telapsed_total: %lf\n", args[i].elapsed_total);
    //     printf("\n\n");
    // }

    if (A_matrix) free(A_matrix);
    if (B_vector) free(B_vector);
    if (X_vector) free(X_vector);
    if (X_indeces) free(X_indeces);
    if (X_true) free(X_true);

    pthread_barrier_destroy(&barrier);
    delete[] args;
    delete[] main_elems;
    free(tids);

    printf("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d p = %d\n\n", argv[0], residual, 
        time_res, formula_type, matrix_size, block_size, threads_amount);

    return 0;
}
