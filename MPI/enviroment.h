#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fenv.h>
#include <math.h>

#include <mpi.h>

// Печать известного участка матрицы (для тестов)
void print_fixed_array(const double *a, int m, int n);
void print_fixed_array(const int *a, int m, int n);

// Получение глобального индекса строки, имея локальный индекс строки в процессе k
int row_local2global(int n, int m, int p, int k, int i_loc);

// Получение локального индекса строки в процессе k, имея глобальный индекс строки
int row_global2local(int n, int m, int p, int k, int i_loc);

// Получения числа блочных строк в матрице (полных и не полных)
int total_blocks(int n, int m);

// Получение максимального числа строк в процессе
int get_max_rows_block(int n, int m, int p);

// Число блочных строк в процессе k
int get_rows_block(int n, int m, int p, int k);

// По глобальному индексу строки получить номер процесса, который ее обслуживает
int get_k(int n, int m, int p, int i_glob);

// По глобальному индексу блочной строки получить локальный индекс блочной строки в процессе k
int block_row_global2local(int n, int m, int p, int k, int block_global);

// Заполнение матрицы по формуле
void init_matrix(double *a, int n, int m, int p, int k, int f_type, double (*f)(int f_type, int size, int i, int j));

// Формула
double formula(int f_type, int size, int i, int j);

// Чтение матрицы из файла
int read_matrix(double *a, int n, int m, int p, int k, const char *name, double *buf /* n*m */);

// Печать матрицы
void print_matrix(const double *a, int n, int m, int p, int k, double *buf, int max_print);

// Печатает из полосы n*m. Не более чем (max_print - printed) - возвращаемое значение.
// max_print - кол-во строк и столбцов
int print_array(const double *a, int n, int m, int printed, int max_print);

// Скопировать блок из буффера
void copyBlockString(const double *String, double *V, int j, int matrix_size, int block_size, int v_size_1, int v_size_2);

// Вставить блок в буффер
void writeBlockString(double *String, const double *V, int j, int matrix_size, int block_size, int v_size_1, int v_size_2);

// Найти норму обратного j-ого блока
double RecipMatrixNorm(double *String, int matrix_size, int block_size, int j, double matrix_A_norm, double *V1, double *V3);

// Найти норму матрицы
double matrixNorm(const double *a, int size_1, int size_2);

// Найти обратную к А матрицу и результат записать в Е
int recipMatrix(double *A, double *E, int m, double matrix_A_norm);

// Сделать матрицу V единичной
void setIdentityMatrix(double *V, int size_1, int size_2);

// Сделать матрицу V нулевой
void setZeroMatrix(double *V, int size_1, int size_2);

// Поменять местами строки матриц А и Е
void changeRows(double *A, double *E, int k, int index, int m);

// Поменять местами блочными колонны с индексами s и main_index в буффере buf и куске матрицы a
void swap_columns(double *a, double *buf, int n, int m, int p, int k, int s, int main_index);

// Перемножить 2 квадратных матрицы размера m
void multiplyMatricesQuadratic(const double *a, const double *b, double *c, int m);

// Перемножить 2 произвольных матрицы
void multiplyMatrices(const double *a, int a_m, int a_n, const double *b, int b_m, int b_n, double *c);

// Скопировать блок матрицы с началом в (i;j)
void copyBlockMatix(double *a, double *V, int i, int j, int n, int m, int v_size_1, int v_size_2);

// Скопировать блок вектора B
void copyBlockVector(const double *b, double *b_buf, int i, int m, int size);

// Записать блок в вектор В
void writeBlockVector(double *b, const double *b_buf, int i, int m, int size);

// Записать блок в матрицу с началом в
void writeBlockMatix(double *a, double *V, int i, int j, int n, int m, int v_size_1, int v_size_2);

// Вычесть из матрицы а матрицу b. Результат записывается в а
void subtractMatrices(double *a, const double *b, int size_1, int size_2);

// Сложить матрицы а и b. Результат записывается в а
void addMatrices(double *a, const double *b, int size_1, int size_2);