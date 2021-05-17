#include "enviroment.h"

int max(int a, int b)
{
    return (a > b ? a : b);
}

int min(int a, int b)
{
    return (a < b ? a : b);
}

void print_fixed_array(const double *a, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.2f ", a[i * n + j]);
        }

        printf("\n");
    }
    printf("\n");
}

void print_fixed_array(const int *a, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%d ", a[i * n + j]);
        }

        printf("\n");
    }
    printf("\n");
}

int row_local2global(int n, int m, int p, int k, int i_loc)
{
    int block_row_loc = i_loc / m;

    (void)n;
    return (k * m + p * block_row_loc * m + i_loc % m);
}

int row_global2local(int n, int m, int p, int k, int i_glob)
{
    int block_row_glob = i_glob / m;
    int block_row_local = block_row_glob / p;

    (void)n;
    (void)k;
    return (block_row_local * m + i_glob % m);
}

int total_blocks(int n, int m)
{
    return (n % m == 0 ? n / m : n / m + 1);
}

int get_max_rows_block(int n, int m, int p)
{
    int blocks = n / m + (n % m == 0 ? 0 : 1); // = (n+(m-1))/m

    return (blocks / p + (blocks % p == 0 ? 0 : 1)); // = (blocks + p-1)/p
}

int get_rows_block(int n, int m, int p, int k)
{
    int blocks = (n + m - 1) / m;

    return ((blocks % p == 0 || blocks % p <= k) ? blocks / p : blocks / p + 1); // (!!исправил!!)
}

int get_k(int n, int m, int p, int i_glob)
{
    (void)n;
    return (i_glob / m) % p;
}

int block_row_global2local(int n, int m, int p, int k, int block_global)
{
    (void)n;
    (void)m;
    (void)k;
    return block_global / p;
}

void init_matrix(double *a, int n, int m, int p, int k, int f_type, double (*f)(int f_type, int size, int i, int j))
{
    /* обычные координаты, не блочные */
    int i_loc, j_loc, i_glob, j_glob, rows_block, rows;

    rows_block = get_rows_block(n, m, p, k); // кол-во блочных строк в процессе k
    rows = rows_block * m;

    for (i_loc = 0; i_loc < rows; i_loc++)
    {
        i_glob = row_local2global(n, m, p, k, i_loc);
        for (j_loc = 0; j_loc < n; j_loc++)
        {
            j_glob = j_loc;

            a[i_loc * n + j_loc] = (*f)(f_type, n, i_glob, j_glob);
        }
    }
}

double formula(int f_type, int size, int i, int j)
{
    double res = 0;

    switch (f_type)
    {
    case 1:
        res = size - max(i, j);
        break;
    case 2:
        res = max(i + 1, j + 1);
        break;
    case 3:
        res = abs(i - j);
        break;
    case 4:
        res = 1.0 / (i + j + 1);
        break;
    default:
        break;
    }

    return res;
}

int read_array(FILE *fp, double *buf, int len, int curr_block, int max_blocks, int last, int n)
{
    int counter = 0;
    int i;
    for (i = 0; i < len; i++)
    {
        if (fscanf(fp, "%lf", buf + i) != 1)
        {
            if (curr_block != max_blocks - 1)
            {
                return -2;
            }
            else if ((last == 1) && (counter % n == 0))
            {
                if (feof(fp) != 0)
                {
                    // файл закончился
                    buf[i] = 0;
                }
                else
                {
                    // какая-то ошибка
                    return -1;
                }
            }
            else
            {
                return -2;
            }
            // if (feof(fp) != 0)
            // {
            //     // файл закончился
            //     for (int j = i; j < len; j++)
            //         buf[j] = 0;
            // }
            // else
            // {
            //     // какая-то ошибка
            //     return 1;
            // }
        }
        else
        {
            counter++;
        }
    }

    return 0;
}

int read_matrix(double *a, int n, int m, int p, int k, const char *name, double *buf /* n*m */)
{
    int main_k = 0;
    FILE *fp = NULL;
    int block;
    int max_blocks = (n + m - 1) / m;
    int err = 0;
    int last = n % m == 0 ? 0 : 1;
    MPI_Comm comm = MPI_COMM_WORLD;

    if (k == main_k)
    {
        fp = fopen(name, "r");
        if (!fp)
            err = 1;
    }
    MPI_Bcast(&err, 1, MPI_INT, main_k, comm);
    if (err > 0)
        return err;

    memset(buf, 0, n * m * sizeof(double));
    for (block = 0; block < max_blocks; block++)
    {
        int owner = block % p;
        int block_loc = block / p;

        // Главный процесс
        if (k == main_k)
        {
            err += read_array(fp, buf, n * m, block, max_blocks, last, n);

            if (owner == main_k)
            {
                // Просто скопировать буффер на место
                memcpy(a + block_loc * m * n, buf, n * m * 8);
                // for (int j = 0; j < m * n; j++)
                //     a[block_loc * m * n + j] = buf[j];
            }
            else
            {
                // Отправить буфер в процесс owner
                MPI_Send(buf, n * m, MPI_DOUBLE, owner, 0 /*tag*/, comm);
            }
        }
        else
        {
            //Остальные процессы
            if (owner == k)
            {
                MPI_Status st;
                MPI_Recv(a + block_loc * m * n, n * m, MPI_DOUBLE, main_k, 0 /*tag*/, comm, &st);
            }
        }
    }

    if (k == main_k)
    {
        fclose(fp);
        fp = 0;
    }

    MPI_Bcast(&err, 1, MPI_INT, main_k, comm);
    if (err != 0)
    {
        printf("error while reading\n");
        return err;
    }

    return 0;
}

void print_matrix(const double *a, int n, int m, int p, int k, double *buf, int max_print)
{
    int main_k = 0;
    int block;
    int max_blocks = (n + m - 1) / m;
    int printed = 0;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status st;
    int block_loc = 0;
    int blocks_in_process = get_rows_block(n, m, p, k);

    for (block = 0; block < max_blocks; block++)
    {
        int owner = block % p;

        if (k == main_k)
        {
            // Главный процесс
            if (owner == main_k)
            {
                // Печатаем свои данные
                printed += print_array(a + block_loc * n * m, n, m, printed, max_print);
                block_loc++;
            }
            else
            {
                MPI_Recv(buf, n * m, MPI_DOUBLE, owner, 0 /*tag*/, comm, &st);
                printed += print_array(buf, n, m, printed, max_print);
            }
        }
        else
        {
            // Остальные процессы
            if (blocks_in_process > block_loc)
            {
                MPI_Send(a + block_loc * n * m, n * m, MPI_DOUBLE, main_k, 0 /*tag*/, comm);
                block_loc++;
            }
        }
    }

    if (k == main_k)
        printf("\n");
}

int print_array(const double *a, int n, int m, int printed, int max_print)
{
    int i, j;
    int can_print = max_print - printed;
    int rows = min(m, can_print);

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < max_print; j++)
        {
            printf("%10.3e ", a[i * n + j]);
        }

        printf("\n");
    }

    return rows;
}

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

void copyBlockVector(const double *b, double *b_buf, int i, int m, int size)
{
    for (int t = 0; t < size; t++)
    {
        b_buf[t] = b[i * m + t];
    }
}

void writeBlockVector(double *b, const double *b_buf, int i, int m, int size)
{
    for (int t = 0; t < size; t++)
    {
        b[i * m + t] = b_buf[t];
    }
}

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

        changeRows(A, E, k, index, m);
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

void changeRows(double *A, double *E, int k, int index, int m)
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

void swap_columns(double *a, double *buf, int n, int m, int p, int k, int s, int main_index)
{
    int i, j;
    double temp;

    int rows_block = get_max_rows_block(n, m, p);
    int rows = rows_block * m;

    // смена строк в куске матрицы а
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < m; j++)
        {
            temp = a[i * n + s * m + j];
            a[i * n + s * m + j] = a[i * n + main_index * m + j];
            a[i * n + main_index * m + j] = temp;
        }
    }

    // смена строк в буффере buf
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < m; j++)
        {
            temp = buf[i * n + s * m + j];
            buf[i * n + s * m + j] = buf[i * n + main_index * m + j];
            buf[i * n + main_index * m + j] = temp;
        }
    }

    (void)k;
}

void copyBlockMatix(double *a, double *V, int i, int j, int n, int m, int v_size_1, int v_size_2)
{
    for (int l = 0; l < v_size_1; l++)
    {
        for (int h = 0; h < v_size_2; h++)
        {
            V[l * v_size_2 + h] = a[i * n * m + j * m + l * n + h];
        }
    }
}

void writeBlockMatix(double *a, double *V, int i, int j, int n, int m, int v_size_1, int v_size_2)
{
    for (int l = 0; l < v_size_1; l++)
    {
        for (int h = 0; h < v_size_2; h++)
        {
            a[i * n * m + j * m + l * n + h] = V[l * v_size_2 + h];
        }
    }
}

void multiplyMatricesQuadratic(const double *a, const double *b, double *c, int m)
{
    double a00, a01, a02,
        a10, a11, a12,
        a20, a21, a22;

    int i, j, k;

    for (i = 0; i < m; i += 3)
    {
        for (j = 0; j < m; j += 3)
        {
            a00 = 0;
            a01 = 0;
            a02 = 0;

            a10 = 0;
            a11 = 0;
            a12 = 0;

            a20 = 0;
            a21 = 0;
            a22 = 0;

            for (k = 0; k < m; k++)
            {
                a00 += a[i * m + k] * b[k * m + j];
                a01 += a[i * m + k] * b[k * m + j + 1];
                a02 += a[i * m + k] * b[k * m + j + 2];

                a10 += a[(i + 1) * m + k] * b[k * m + j];
                a11 += a[(i + 1) * m + k] * b[k * m + j + 1];
                a12 += a[(i + 1) * m + k] * b[k * m + j + 2];

                a20 += a[(i + 2) * m + k] * b[k * m + j];
                a21 += a[(i + 2) * m + k] * b[k * m + j + 1];
                a22 += a[(i + 2) * m + k] * b[k * m + j + 2];
            }

            c[i * m + j] = a00;
            c[i * m + j + 1] = a01;
            c[i * m + j + 2] = a02;

            c[(i + 1) * m + j] = a10;
            c[(i + 1) * m + j + 1] = a11;
            c[(i + 1) * m + j + 2] = a12;

            c[(i + 2) * m + j] = a20;
            c[(i + 2) * m + j + 1] = a21;
            c[(i + 2) * m + j + 2] = a22;
        }
    }
}

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