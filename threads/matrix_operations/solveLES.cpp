#include "matrix_operations.h"

///////////////////////////////  начало ОБЪЯВЛЕНИЯ INLINE-ФУНКЦИЙ  ///////////////////////////////

void multiplyMatricesQuadratic(const double *a, const double *b, double *c, int m);
void copyBlockMatrix(const double *A, double *V, int i, int j, int matrix_size, int block_size, int v_size_1, int v_size_2);
void copyBlockString(const double *String, double *V, int j, int matrix_size, int block_size, int v_size_1, int v_size_2);
void writeBlockMatrix(double *A, const double *V, int i, int j, int matrix_size, int block_size, int v_size_1, int v_size_2);
void subtractMatrices(double *a, const double *b, int size_1, int size_2);

inline void multiplyMatricesQuadratic(const double *a, const double *b, double *c, int m) {
  double a00, a01, a02,
         a10, a11, a12,
         a20, a21, a22;
  
  int i, j, k;

  for (i=0; i<m; i+=3) {
    for (j=0; j<m; j+=3) {
      a00 = 0;
      a01 = 0;
      a02 = 0;

      a10 = 0;
      a11 = 0;
      a12 = 0;

      a20 = 0;
      a21 = 0;
      a22 = 0;

      for (k=0; k<m; k++) {
        a00 += a[i*m + k] * b[k*m + j];
        a01 += a[i*m + k] * b[k*m + j+1];
        a02 += a[i*m + k] * b[k*m + j+2];

        a10 += a[(i+1)*m + k] * b[k*m + j];
        a11 += a[(i+1)*m + k] * b[k*m + j+1];
        a12 += a[(i+1)*m + k] * b[k*m + j+2];

        a20 += a[(i+2)*m + k] * b[k*m + j];
        a21 += a[(i+2)*m + k] * b[k*m + j+1];
        a22 += a[(i+2)*m + k] * b[k*m + j+2];
      }

      c[i*m + j] = a00;
      c[i*m + j+1] = a01;
      c[i*m + j+2] = a02;

      c[(i+1)*m + j] = a10;
      c[(i+1)*m + j+1] = a11;
      c[(i+1)*m + j+2] = a12;

      c[(i+2)*m + j] = a20;
      c[(i+2)*m + j+1] = a21;
      c[(i+2)*m + j+2] = a22;
    }
  }
}

inline void copyBlockMatrix(const double *A, double *V, int i, int j, int matrix_size, int block_size, int v_size_1, int v_size_2) {
  /*
    копирование из блока с адресом (i;j) 
    v_size_1 и v_size_2 - стандартно равны block_size, но могут так же принимать значение l - остатка
  */
  for (int t=0; t<v_size_1; t++) {
    for (int s=0; s<v_size_2; s++) {
      V[t*v_size_2 + s] = A[i*matrix_size*block_size + j*block_size + t*matrix_size + s];
    }
  }
}

inline void copyBlockString(const double *String, double *V, int j, int matrix_size, int block_size, int v_size_1, int v_size_2) {
  for (int t=0; t<v_size_1; t++) {
    for (int k=0; k<v_size_2; k++) {
      V[t*v_size_2 + k] = String[t*matrix_size + j*block_size + k];
    }
  }
}

inline void writeBlockMatrix(double *A, const double *V, int i, int j, int matrix_size, int block_size, int v_size_1, int v_size_2) {
  /* 
    запись в блок с адресом (i;j)
    v_size_1 и v_size_2 - стандартно равны block_size, но могут так же принимать значение l - остатка  
  */
  for (int t=0; t<v_size_1; t++) {
    for (int s=0; s<v_size_2; s++) {
      A[i*matrix_size*block_size + j*block_size + t*matrix_size + s] = V[t*v_size_2 + s];
    }
  }
}

inline void subtractMatrices(double *a, const double *b, int size_1, int size_2) {  
  for (int i=0; i<size_1; i++) {
    for (int j=0; j<size_2; j++) {
      a[i*size_2 + j] = a[i*size_2 + j] - b[i*size_2 + j];
    }
  }
}



///////////////////////////////  КОНЕЦ ОБЪЯВЛЕНИЯ INLINE-ФУНКЦИЙ  ///////////////////////////////

int solveLES(double *A, int matrix_size, int block_size, double *B, int *X_indeces, pthread_barrier_t *barrier,
    int t_amount, int t_index, int formula_type, mainElem *main_elems)
    {
    double *V1 = NULL;
    double *V2 = NULL;
    double *V3 = NULL;

    double *B1 = NULL;
    double *B2 = NULL;
    double *B3 = NULL;

    double *String = NULL;

    // различные счетчики
    int k = 0;
    int i = 0;
    int h = 0;
    int j = 0;

    // размер остаточного блока
    int l = 0;

    // индикатор того, что размер матрицы не поделился на размер блока
    int last = 0;

    // индикатор того, что обратная матрица найдена
    int success = 0;

    int p = t_amount;

    // индекс в строке или индекс строки, с которого поток начинает свою работу
    int start_point = 0;
    int main_index = 0;
    int main_found = -1;
    int swap_index = 0;
    int recip_found = -1;
    int temp;
    int method_works = 1;

    double matrix_A_norm = 0;
    double loc_res = 0;
    double min_norm = 0;

    k = matrix_size/block_size;
    l =  matrix_size - block_size*k;

    V1 = (double*)calloc(block_size*block_size, sizeof(double));
    V2 = (double*)calloc(block_size*block_size, sizeof(double));
    V3 = (double*)calloc(block_size*block_size, sizeof(double));

    B1 = (double*)calloc(block_size, sizeof(double));
    B2 = (double*)calloc(block_size, sizeof(double));
    B3 = (double*)calloc(block_size, sizeof(double));

    String = (double*)calloc(block_size*matrix_size, sizeof(double));

    last = (l>0 ? 1 : 0);

    if (formula_type != 0) {
        for (h=t_index; h<k; h+=p) {
            for (i=0; i<block_size; i++) {
                for (j=0; j<matrix_size; j++) {
                    A[h*block_size*matrix_size + i*matrix_size + j] = 0;
                }
            }
        }
        if (last > 0 && t_index == 0) {
            for (i=0; i<l; i++) {
                for (j=0; j<matrix_size; j++) {
                    A[k*block_size*matrix_size + i*matrix_size + j] = 0;
                }
            }
        }

        /* Проинициализировали все переменные */

        for (h=t_index; h<k; h+=p) {
            for (i=0; i<block_size; i++) {
                for (j=0; j<matrix_size; j++) {
                    A[h*block_size*matrix_size + i*matrix_size + j] = buildMatrix(formula_type, matrix_size, h*block_size + i, j);
                }

                B[h*block_size + i] = sum(A, matrix_size, h*block_size + i);
            }
        }
        if (last > 0 && t_index == 0) {
            for (i=0; i<l; i++) {
                for (j=0; j<matrix_size; j++) {
                    A[k*block_size*matrix_size + i*matrix_size + j] = buildMatrix(formula_type, matrix_size, k*block_size + i, j);
                }

                B[k*block_size + i] = sum(A, matrix_size, k*block_size + i);
            }
        }

        pthread_barrier_wait(barrier);
    }

    // if (t_index == 0) printMatrix(A, matrix_size, matrix_size, matrix_size);
    pthread_barrier_wait(barrier);
    // if (t_index == 0) printMatrix(B, matrix_size, 1, matrix_size);

    matrix_A_norm = matrixNorm(A, matrix_size, matrix_size);


    /*
     ________________________________________________ПРЯМОЙ ХОД МЕТОДА ГАУССА________________________________________________ 
    */

    for (i=0; i<k; i++) {
        /* Внесем в память узла строку нашего шага */
        for(int t=0; t<block_size; t++) {
            for (j = i*block_size; j<matrix_size; j++) {
                String[t*matrix_size + j] = A[i*block_size*matrix_size + t*matrix_size + j];
            }
        }
        // if (t_index == p/2) printMatrix(String, block_size, matrix_size, matrix_size);

        min_norm = BIG_NUMBER;
        main_found = -1;
        swap_index = i;
        recip_found = -1;
        
        /* Перед каждым шагом метода гаусса нужно выбрать гланвый элемент по строке */

        start_point = t_index;
        while (start_point < i) {
            start_point += p;
        }

        for (j=start_point; j<k; j+=p) {
            /* вернет -1 в случае всех вырожденного в строке,норму в случае успеха */
            loc_res = RecipMatrixNorm(String, matrix_size, block_size, j, matrix_A_norm, V1, V3);
            if (loc_res > 0) {
                if (min_norm > loc_res) {
                    min_norm = loc_res;
                    main_index = j;
                }

                main_found = 1;
            }
        }
        if (main_found > 0) {
            main_elems[t_index].index = main_index;
            main_elems[t_index].value = min_norm;
            main_elems[t_index].error = 1;
        }
        else {
            main_elems[t_index].error = -1;
        }

        if (p == 1) {
            if (main_found > 0) {
                recip_found = 1;
                swap_index = main_elems[t_index].index;
            }
        }

        pthread_barrier_wait(barrier);

        reduce_sum(p, main_elems, &swap_index, &recip_found);

        if (recip_found < 0) {
            // Алгоритм нерпименим, запишем единичные матрицы чтобы сэмитировать работу программы
            setIdentityMatrix(V1, block_size, block_size);

            writeBlockString(String, V1, i, matrix_size, block_size, block_size, block_size);

            if (t_index == 0) writeBlockMatrix(A, V1, i, i, matrix_size, block_size, block_size, block_size); 

            method_works = -1;
        } 
        else if ( swap_index > i ) {
            // printf("swap: %d and %d\n", i, swap_index);
            for (h=t_index; h<k; h+=p) {
                changeBlocks(A, matrix_size, block_size, V1, V2, block_size, block_size, h, i, h, swap_index);
            }

            if (last > 0 && t_index == 0 ) {
                changeBlocks(A, matrix_size, block_size, V1, V2, l, block_size, k, i, k, swap_index);
            }

            swapBlockString(String, V1, V2, swap_index, i, matrix_size, block_size);

            // Поменять индексы
            if (t_index == 0) {
                for (int t=0; t<block_size; t++) {
                    temp = X_indeces[i*block_size + t];

                    X_indeces[i*block_size + t] = X_indeces[swap_index*block_size + t];
                    X_indeces[swap_index*block_size + t] = temp;
                }
            }
        }

        pthread_barrier_wait(barrier);

        copyBlockString(String, V1, i, matrix_size, block_size, block_size, block_size);
        
        recipMatrix(V1, V2, block_size, matrix_A_norm);

        /* обратная лежит в V2 */
        if (start_point == i) {
            setIdentityMatrix(V1, block_size, block_size);
            writeBlockMatrix(A, V1, i, i, matrix_size, block_size, block_size, block_size);
            start_point += p;
        }

        for (h=start_point; h<k; h+=p) {
            copyBlockString(String, V1, h, matrix_size, block_size, block_size, block_size);
            multiplyMatricesQuadratic(V2, V1, V3, block_size);

            writeBlockMatrix(A, V3, i, h, matrix_size, block_size, block_size, block_size);
        }
        if (last > 0 && t_index == 0) {
            copyBlockString(String, V1, k, matrix_size, block_size, block_size, l);
            multiplyMatrices(V2, block_size, block_size, V1, block_size, l, V3);

            writeBlockMatrix(A, V3, i, k, matrix_size, block_size, block_size, l);
        }
        if (t_index == 0) {
            copyBlockVector(B, B1, i, block_size, block_size);
            multiplyMatrices(V2, block_size, block_size, B1, block_size, 1, B2);
            writeBlockVector(B, B2, i, block_size, block_size);
        }

        pthread_barrier_wait(barrier);

        /* Обновляем строку */
        for(int t=0; t<block_size; t++) {
            for (j = i*block_size; j<matrix_size; j++) {
                String[t*matrix_size + j] = A[i*block_size*matrix_size + t*matrix_size + j];
            }
        }

        // if (t_index == p/2) printMatrix(A, matrix_size, matrix_size, matrix_size);

        /* ВЫЧИТАНИЕ СТРОК: каждый поток вычитает из своего набора строк */

        for (int t=start_point; t<k; t+=p) {
            // if (t == i) continue;
            copyBlockMatrix(A, V1, t, i, matrix_size, block_size, block_size, block_size);

            for (j=i; j<k; j++) {
                copyBlockString(String, V2, j, matrix_size, block_size, block_size, block_size);
                multiplyMatricesQuadratic(V1, V2, V3, block_size);

                copyBlockMatrix(A, V2, t, j, matrix_size, block_size, block_size, block_size);
                subtractMatrices(V2, V3, block_size, block_size);

                writeBlockMatrix(A, V2, t, j, matrix_size, block_size, block_size, block_size);
            }

            if (last > 0) {
                copyBlockString(String, V2, k, matrix_size, block_size, block_size, l);
                multiplyMatrices(V1, block_size, block_size, V2, block_size, l, V3);

                copyBlockMatrix(A, V2, t, k, matrix_size, block_size, block_size, l);
                subtractMatrices(V2, V3, block_size, l);

                writeBlockMatrix(A, V2, t, k, matrix_size, block_size, block_size, l);
            }

            copyBlockVector(B, B1, i, block_size, block_size);
            multiplyMatrices(V1, block_size, block_size, B1, block_size, 1, B2);

            copyBlockVector(B, B1, t, block_size, block_size);
            subtractMatrices(B1, B2, block_size, 1);

            writeBlockVector(B, B1, t, block_size, block_size);
        }

        if (last > 0 && t_index == 0) {
            copyBlockMatrix(A, V1, k, i, matrix_size, block_size, l, block_size);

            for (j=i; j<k; j++) {
                copyBlockString(String, V2, j, matrix_size, block_size, block_size, block_size);
                multiplyMatrices(V1, l, block_size, V2, block_size, block_size, V3);

                copyBlockMatrix(A, V2, k, j, matrix_size, block_size, l, block_size);
                subtractMatrices(V2, V3, l, block_size);

                writeBlockMatrix(A, V2, k, j, matrix_size, block_size, l, block_size);
            }

            copyBlockString(String, V2, k, matrix_size, block_size, block_size, l);
            multiplyMatrices(V1, l, block_size, V2, block_size, l, V3);

            copyBlockMatrix(A, V2, k, k, matrix_size, block_size, l, l);
            subtractMatrices(V2, V3, l, l);

            writeBlockMatrix(A, V2, k, k, matrix_size, block_size, l, l);

            /* Нижний блок веткора В */
            copyBlockVector(B, B1, i, block_size, block_size);
            multiplyMatrices(V1, l, block_size, B1, block_size, 1, B2);

            copyBlockVector(B, B1, k, block_size, l);
            subtractMatrices(B1, B2, l, 1);

            writeBlockVector(B, B1, k, block_size, l);
        }

        pthread_barrier_wait(barrier);

    }

    // Обработка остаточного элемента в правом нижнем углу
    if (last > 0 && t_index == 0) {
        copyBlockMatrix(A, V1, k, k, matrix_size, block_size, l, l);
        success = recipMatrix(V1, V2, l, matrix_A_norm);

        if ( !success ) {
            // Алгоритм нерпименим, запишем единичные матрицы чтобы сэмитировать работу программы
            setIdentityMatrix(V2, l, l);

            writeBlockMatrix(A, V2, k, k, matrix_size, block_size, l, l);

            method_works = -1;
        }

        copyBlockVector(B, B1, k, block_size, l);
        multiplyMatrices(V2, l, l, B1, l, 1, B2);

        writeBlockVector(B, B2, k, block_size, l);

        setIdentityMatrix(V1, l, l);
        writeBlockMatrix(A, V1, k, k, matrix_size, block_size, l, l);
    }

    pthread_barrier_wait(barrier);

    //  ________________________________________________ОБРАТНЫЙ ХОД МЕТОДА ГАУССА________________________________________________ 
    start_point = t_index;
    
    // Обработка последнего столбца
    if (last > 0) {
        copyBlockVector(B, B1, k, block_size, l);

        for (int t=start_point; t<k; t+=p) {
            copyBlockMatrix(A, V1, t, k, matrix_size, block_size, block_size, l);
            multiplyMatrices(V1, block_size, l, B1, l, 1, B2);

            copyBlockVector(B, B3, t, block_size, block_size);
            subtractMatrices(B3, B2, block_size, 1);

            writeBlockVector(B, B3, t, block_size, block_size);

        }
    }

    pthread_barrier_wait(barrier);

    for (i=1; i<k; i++) {
        copyBlockVector(B, B1, k-i, block_size, block_size);

        if (start_point == (k-i)) start_point+=p; 

        for (int t=start_point; t<k; t+=p) {
            if (t == (k-i)) continue;
            
            copyBlockMatrix(A, V1, t, k-i, matrix_size, block_size, block_size, block_size);
            multiplyMatrices(V1, block_size, block_size, B1, block_size, 1, B2);

            copyBlockVector(B, B3, t, block_size, block_size);
            subtractMatrices(B3, B2, block_size, 1);

            writeBlockVector(B, B3, t, block_size, block_size);

        }

        pthread_barrier_wait(barrier);
    }

    if (V1) free(V1);
    if (V2) free(V2);
    if (V3) free(V3);
    if (B1) free(B1);
    if (B2) free(B2);
    if (B3) free(B3);  

    if (String) free(String);

    if (method_works > 0) return 1;
    else return -1;
}