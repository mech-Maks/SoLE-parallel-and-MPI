#include "matrix_operations.h"

/* For build */
int max(int a, int b)
{
  return (a > b ? a : b);
}

int module(int a)
{
  return (a < 0 ? -a : a);
}

double sum(const double *A, int size, int index)
{
  double res;
  int limit = (size + 1) / 2;
  res = 0;

  for (int j = 0; j < limit; j++)
  {
    res += A[index * size + 2 * j];
  }

  return res;
}

int buildMatrixFromFile(double *A, int size, int type, const char *filename, double *B)
{
  FILE *from = NULL;
  double c = 0;
  int i = 0;
  int temp = 0;
  int status = 0;

  switch (type)
  {
  case 0:
    i = 0;
    temp = size * size;
    from = fopen(filename, "r");
    if (from == NULL)
    {
      printf("Can't open %s\n", filename);
      status = 1;
      break;
    }

    while ((fscanf(from, "%lf", &c) == 1) && (i < temp))
    {
      A[i] = c;
      i++;
    }

    if (i != temp)
    {
      printf("Not numbers in file or not enough numbers\n");
      status = 1;
    }
    fclose(from);
    break;

  default:
    break;
    status = 1;
  }

  if (!status)
  {
    for (i = 0; i < size; i++)
    {
      B[i] = sum(A, size, i);
    }

    return 1;
  }

  return 0;
}

double buildMatrix(int type, int size, int i, int j)
{
  double res = 0;

  switch (type)
  {
  case 1:
    res = size - max(i, j);
    break;
  case 2:
    res = max(i + 1, j + 1);
    break;
  case 3:
    res = module(i - j);
    break;
  case 4:
    res = 1.0 / (i + j + 1);
    break;
  default:
    break;
  }

  return res;
}

/* For inverse */
int recipMatrix(double *A, double *E, int m, double matrix_A_norm)
{
  int i = 0;
  int j = 0;
  int k = 0;
  int index = 0;

  double main_elem = 0;
  double temp = 0;

  setIdentityMatrix(E, m, m);
  while (k < m)
  {
    main_elem = fabs(A[k * m + k]);
    index = k;

    for (i = k + 1; i < m; i++)
    {
      temp = fabs(A[i * m + k]);

      if (temp > main_elem)
      {
        main_elem = temp;
        index = i;
      }
    }

    changeRaws(A, E, k, index, m);
    main_elem = A[k * m + k];

    if (fabs(main_elem) < 1e-16 * matrix_A_norm)
    {
      return 0;
    }

    temp = 1.0 / main_elem;
    for (j = 0; j < m; j++)
    {
      A[k * m + j] = (double)A[k * m + j] * temp;
      E[k * m + j] = (double)E[k * m + j] * temp;
    }
    // A[k*m + k] = 1;

    for (i = k + 1; i < m; i++)
    {
      temp = A[i * m + k];
      for (j = 0; j < m; j++)
      {
        A[i * m + j] -= A[k * m + j] * temp;
        E[i * m + j] -= E[k * m + j] * temp;
      }
      // A[i*m + k] = 0;
    }

    k++;
  }

  k--;
  temp = A[k * m + k];
  for (j = 0; j < m; j++)
  {
    A[k * m + j] = (double)A[k * m + j] / temp;
    E[k * m + j] = (double)E[k * m + j] / temp;
  }

  while (k > 0)
  {
    for (i = k - 1; i >= 0; i--)
    {
      temp = A[i * m + k];
      for (j = 0; j < m; j++)
      {
        A[i * m + j] = A[i * m + j] - A[k * m + j] * temp;
        E[i * m + j] = E[i * m + j] - E[k * m + j] * temp;
      }
    }

    k--;
  }

  return 1;
}

void changeRaws(double *A, double *E, int k, int index, int m)
{
  double temp = 0;

  for (int i = 0; i < m; i++)
  {
    temp = A[index * m + i];
    A[index * m + i] = A[k * m + i];
    A[k * m + i] = temp;

    temp = E[index * m + i];
    E[index * m + i] = E[k * m + i];
    E[k * m + i] = temp;
  }
}

void setIdentityMatrix(double *V, int size_1, int size_2)
{
  for (int i = 0; i < size_1; i++)
  {
    for (int j = 0; j < size_2; j++)
    {
      if (i != j)
      {
        V[i * size_2 + j] = 0;
      }
      else
      {
        V[i * size_2 + j] = 1;
      }
    }
  }
}

void setZeroMatrix(double *V, int size_1, int size_2)
{
  for (int i = 0; i < size_1; i++)
  {
    for (int j = 0; j < size_2; j++)
    {
      V[i * size_2 + j] = 0;
    }
  }
}

/* For norm */
double matrixNorm(const double *a, int size_1, int size_2)
{
  double max = 0;
  double temp = 0;

  for (int j = 0; j < size_2; j++)
  {
    temp = 0;

    for (int i = 0; i < size_1; i++)
    {
      temp += fabs(a[i * size_2 + j]);
    }

    if (temp > max)
    {
      max = temp;
    }
  }

  return max;
}

double discrepancy(const double *A, const double *x, const double *b, int size)
{
  double *c = NULL;

  double norm1 = 0;
  double norm2 = 0;
  double res = 0;

  c = (double *)malloc(size * sizeof(double));
  setZeroMatrix(c, size, 1);

  // printMatrix(A, size, size, size);
  // printMatrix(x, size, 1, size);

  multiplyMatrices(
      A, size, size,
      x, size, 1,
      c);

  // printf("\tc:\n");
  // printMatrix(c, size, 1, size);

  for (int i = 0; i < size; i++)
  {
    // printf("c[%d]: %e\n", i, c[i]);
    // printf("b[%d]: %e\n", i, b[i]);
    c[i] -= b[i];
  }

  norm1 = matrixNorm(c, size, 1);
  norm2 = matrixNorm(b, size, 1);

  res = (double)1 / norm2;

  // printf("%e\n", res);
  // printf("%e\n", norm1);

  res = res * norm1;

  if (c)
    free(c);

  return res;
}

/* For multiplication */
void multiplyMatrices(
    const double *a, int a_m, int a_n,
    const double *b, int b_m, int b_n,
    double *c)
{
  // if (a_n != b_m) {
  //     /* Проверка для себя */
  //     printf("Can't multiply. Bad matrices\n");
  //     exit(-2);
  // }

  int i = 0;
  int j = 0;
  int k = 0;

  double temp = 0;
  double temp1 = 0;
  double temp2 = 0;
  double temp3 = 0;
  double temp4 = 0;

  int last = a_m % 4;

  /* если кол-во строк в матрице А кратно 4 - все отлично */
  for (i = 0; i < a_m - 3; i += 4)
  {
    for (j = 0; j < b_n; j++)
    {
      temp1 = 0;
      temp2 = 0;
      temp3 = 0;
      temp4 = 0;

      for (k = 0; k < b_m; k++)
      {
        temp = b[k * b_n + j];
        temp1 = temp1 + a[i * a_n + k] * temp;
        temp2 = temp2 + a[(i + 1) * a_n + k] * temp;
        temp3 = temp3 + a[(i + 2) * a_n + k] * temp;
        temp4 = temp4 + a[(i + 3) * a_n + k] * temp;
      }

      c[i * b_n + j] = temp1;
      c[(i + 1) * b_n + j] = temp2;
      c[(i + 2) * b_n + j] = temp3;
      c[(i + 3) * b_n + j] = temp4;
    }
  }

  /* краевой случай, если в матрице А кол-во строк не делится на 4 */
  if (last != 0)
  {
    i = a_m - last;

    while (last > 0)
    {
      for (j = 0; j < b_n; j++)
      {
        temp1 = 0;

        for (k = 0; k < b_m; k++)
        {
          temp1 = temp1 + a[(i + last - 1) * a_n + k] * b[k * b_n + j];
        }

        c[(i + last - 1) * b_n + j] = temp1;
      }

      last--;
    }
  }
}

// void multiplyMatricesQuadratic(const double *a, const double *b, double *c, int m) {
//   double a00, a01, a02,
//          a10, a11, a12,
//          a20, a21, a22;

//   int i, j, k;

//   for (i=0; i<m; i+=3) {
//     for (j=0; j<m; j+=3) {
//       a00 = 0;
//       a01 = 0;
//       a02 = 0;

//       a10 = 0;
//       a11 = 0;
//       a12 = 0;

//       a20 = 0;
//       a21 = 0;
//       a22 = 0;

//       for (k=0; k<m; k++) {
//         a00 += a[i*m + k] * b[k*m + j];
//         a01 += a[i*m + k] * b[k*m + j+1];
//         a02 += a[i*m + k] * b[k*m + j+2];

//         a10 += a[(i+1)*m + k] * b[k*m + j];
//         a11 += a[(i+1)*m + k] * b[k*m + j+1];
//         a12 += a[(i+1)*m + k] * b[k*m + j+2];

//         a20 += a[(i+2)*m + k] * b[k*m + j];
//         a21 += a[(i+2)*m + k] * b[k*m + j+1];
//         a22 += a[(i+2)*m + k] * b[k*m + j+2];
//       }

//       c[i*m + j] = a00;
//       c[i*m + j+1] = a01;
//       c[i*m + j+2] = a02;

//       c[(i+1)*m + j] = a10;
//       c[(i+1)*m + j+1] = a11;
//       c[(i+1)*m + j+2] = a12;

//       c[(i+2)*m + j] = a20;
//       c[(i+2)*m + j+1] = a21;
//       c[(i+2)*m + j+2] = a22;
//     }
//   }
// }

/* For matrix sum */
void addMatrices(double *a, const double *b, int size_1, int size_2)
{
  for (int i = 0; i < size_1; i++)
  {
    for (int j = 0; j < size_2; j++)
    {
      a[i * size_2 + j] = a[i * size_2 + j] + b[i * size_2 + j];
    }
  }
}

/* result will be in massive a */
void subtractMatrices(double *a, const double *b, int size_1, int size_2)
{
  for (int i = 0; i < size_1; i++)
  {
    for (int j = 0; j < size_2; j++)
    {
      a[i * size_2 + j] = a[i * size_2 + j] - b[i * size_2 + j];
    }
  }
}

/* For printing */
void printMatrix(const double *A, int size_1, int size_2, int max_size)
{
  int m = size_1;
  int n = size_2;

  if (m > max_size)
  {
    m = max_size;
  }

  if (n > max_size)
  {
    n = max_size;
  }

  printf("\n");
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf(" %10.3e", A[i * size_2 + j]);
    }
    printf("\n\n");
  }
  printf("\n");
}

/* Get block from matrix */
void copyBlockMatrix(const double *A, double *V, int i, int j, int matrix_size, int block_size, int v_size_1, int v_size_2)
{
  /*
    копирование из блока с адресом (i;j) 
    v_size_1 и v_size_2 - стандартно равны block_size, но могут так же принимать значение l - остатка
  */
  for (int t = 0; t < v_size_1; t++)
  {
    for (int s = 0; s < v_size_2; s++)
    {
      V[t * v_size_2 + s] = A[i * matrix_size * block_size + j * block_size + t * matrix_size + s];
    }
  }
}

/* Write in block */
void writeBlockMatrix(double *A, const double *V, int i, int j, int matrix_size, int block_size, int v_size_1, int v_size_2)
{
  /* 
    запись в блок с адресом (i;j)
    v_size_1 и v_size_2 - стандартно равны block_size, но могут так же принимать значение l - остатка  
  */
  for (int t = 0; t < v_size_1; t++)
  {
    for (int s = 0; s < v_size_2; s++)
    {
      A[i * matrix_size * block_size + j * block_size + t * matrix_size + s] = V[t * v_size_2 + s];
    }
  }
}

void changeBlocks(double *A, int matrix_size, int block_size, double *V1, double *V2, int v_size_1, int v_size_2,
                  int change_pos_1_i, int change_pos_1_j, int change_pos_2_i, int change_pos_2_j)
{
  /* Скопировали блок размера v_size_1 x v_size_2 из блока по адресу (change_pos_1_i; change_pos_1_j) */
  copyBlockMatrix(A, V1, change_pos_1_i, change_pos_1_j, matrix_size, block_size, v_size_1, v_size_2);

  /* Скопировали блок размера v_size_1 x v_size_2 из блока по адресу (change_pos_2_i; change_pos_2_j) */
  copyBlockMatrix(A, V2, change_pos_2_i, change_pos_2_j, matrix_size, block_size, v_size_1, v_size_2);

  /* 
      Записываем в блок размера v_size_1 x v_size_2 по адресу (change_pos_1_i; change_pos_1_j) 
      блок того же размера с адресом (change_pos_2_i; change_pos_2_j)
    */
  writeBlockMatrix(A, V2, change_pos_1_i, change_pos_1_j, matrix_size, block_size, v_size_1, v_size_2);

  /* 
      Записываем в блок размера v_size_1 x v_size_2 по адресу (change_pos_2_i; change_pos_2_j) 
      блок того же размера с адресом (change_pos_1_i; change_pos_1_j)
    */
  writeBlockMatrix(A, V1, change_pos_2_i, change_pos_2_j, matrix_size, block_size, v_size_1, v_size_2);
}

/* Copy block from vector */
void copyBlockVector(const double *B, double *V, int i, int block_size, int v_size)
{
  for (int j = 0; j < v_size; j++)
  {
    V[j] = B[i * block_size + j];
  }
}

/* write block in vector */
void writeBlockVector(double *B, const double *V, int i, int block_size, int v_size)
{
  for (int j = 0; j < v_size; j++)
  {
    B[i * block_size + j] = V[j];
  }
}

/* 
  Find main block in i row and swap it with current diagonal block in i-step 
  1 - success, 0 - fault
*/
double RecipMatrixNorm(double *String, int matrix_size, int block_size, int j, double matrix_A_norm, double *V1, double *V3)
{
  double temp = 0;

  int success = 0;

  copyBlockString(String, V1, j, matrix_size, block_size, block_size, block_size);
  success = recipMatrix(V1, V3, block_size, matrix_A_norm);

  if (success == 1)
  {
    temp = matrixNorm(V3, block_size, block_size);
    return temp;
  }
  else
  {
    return -1;
  }

  return -1;
}

void reduce_sum(int p, mainElem *main_elems, int *swap_index, int *recip_found)
{
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;

  static mainElem *elems = 0;

  static int t_in = 0;
  static int t_out = 0;

  // локальные переменные
  double temp = INFINITY;

  if (p <= 1)
    return;

  pthread_mutex_lock(&mutex);

  if (!elems)
    elems = main_elems;
  else
    for (int i = 0; i < p; i++)
    {
      elems[i].index = main_elems[i].index;
      elems[i].value = main_elems[i].value;
      elems[i].error = main_elems[i].error;
    }

  t_in++;
  if (t_in >= p)
  {
    for (int i = 0; i < p; i++)
    {
      elems[i].index = main_elems[i].index;
      elems[i].value = main_elems[i].value;
      elems[i].error = main_elems[i].error;
    }
    t_out = 0;
    pthread_cond_broadcast(&c_in);
  }
  else
    while (t_in < p)
      pthread_cond_wait(&c_in, &mutex);
  // все вошли

  for (int i = 0; i < p; i++)
  {
    main_elems[i].index = elems[i].index;
    main_elems[i].value = elems[i].value;
    main_elems[i].error = elems[i].error;
  }

  for (int i = 0; i < p; i++)
  {
    if (elems[i].error < 0)
      continue;
    else
    {
      if (elems[i].value < temp)
      {
        temp = elems[i].value;

        (*recip_found) = 1;
        (*swap_index) = i;
      }
    }
  }

  t_out++;
  if (t_out >= p)
  {
    t_in = 0;
    pthread_cond_broadcast(&c_out);
  }
  else
    while (t_out < p)
      pthread_cond_wait(&c_out, &mutex);
  // все вышли

  pthread_mutex_unlock(&mutex);
}

void *thread_func(void *arg)
{
  arguments *args = (arguments *)arg;
  int t_index = args->t_index;

  cpu_set_t cpu;
  int nprocs = get_nprocs();
  CPU_ZERO(&cpu);
  CPU_SET(nprocs - 1 - t_index, &cpu);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu), &cpu);

  pthread_barrier_t *barrier = args->barrier;

  pthread_barrier_wait(barrier);
  // создали все потоки
  args->elapsed_thread = get_time();
  args->elapsed_total = get_full_time();

  /* 1 в случае успеха, -1 в случае краша */
  args->method_does_not_work = solveLES(
      args->A, args->matrix_size, args->block_size, args->B, args->X_indeces, barrier,
      args->t_amount, args->t_index, args->formula_type, args->main_elems);
  args->elapsed_thread = get_time() - args->elapsed_thread;
  args->elapsed_total = get_full_time() - args->elapsed_total;
  pthread_barrier_wait(barrier);
  return 0;
}

/* work with time */
double get_time()
{
  struct rusage r;
  getrusage(RUSAGE_THREAD, &r);
  return r.ru_utime.tv_sec + r.ru_utime.tv_usec / 1.e6;
}

double get_full_time()
{
  struct timeval r;
  gettimeofday(&r, 0);
  return r.tv_sec + r.tv_usec / 1.e6;
}

/* Copy block from String */
void copyBlockString(const double *String, double *V, int j, int matrix_size, int block_size, int v_size_1, int v_size_2)
{
  for (int t = 0; t < v_size_1; t++)
  {
    for (int k = 0; k < v_size_2; k++)
    {
      V[t * v_size_2 + k] = String[t * matrix_size + j * block_size + k];
    }
  }
}

/* write block to string */
void writeBlockString(double *String, const double *V, int j, int matrix_size, int block_size, int v_size_1, int v_size_2)
{
  for (int t = 0; t < v_size_1; t++)
  {
    for (int k = 0; k < v_size_2; k++)
    {
      String[t * matrix_size + j * block_size + k] = V[t * v_size_2 + k];
    }
  }
}

void swapBlockString(double *String, double *V1, double *V2, int index_1, int index_2, int matrix_size, int block_size)
{
  copyBlockString(String, V1, index_1, matrix_size, block_size, block_size, block_size);
  copyBlockString(String, V2, index_2, matrix_size, block_size, block_size, block_size);

  writeBlockString(String, V1, index_2, matrix_size, block_size, block_size, block_size);
  writeBlockString(String, V2, index_1, matrix_size, block_size, block_size, block_size);
}