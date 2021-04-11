#pragma once

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <time.h>
#include <fenv.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <sched.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>

#define BIG_NUMBER 1.e10

class mainElem {
  public:
    double value;

    int index;
    int error;
  
  mainElem() {
    error = 1;
    index = 0;
    value = 0.0;
  }

  
};

class arguments {
  public:
    double *A;
    double *B;

    double elapsed_thread;
    double elapsed_total;

    int *X_indeces;

    int matrix_size;
    int block_size;
    int formula_type;
    int t_index;
    int t_amount;

    int method_does_not_work;

    pthread_barrier_t *barrier;

    mainElem *main_elems;

  arguments() = default;

  
};

/* For build */
int max(int a, int b);
int module(int a);
double sum (const double *A, int size, int index);
int buildMatrixFromFile(double *A, int size, int type, const char *filename, double *B);
double buildMatrix(int type, int size, int i, int j);

/* 
  For inverse 
  1 - success, 0 - fault
*/
int recipMatrix(double *A, double *E, int m, double matrix_A_norm);
void changeRaws(double *A, double *E, int k, int index, int m);
void setIdentityMatrix(double *V, int size_1, int size_2);
void setZeroMatrix(double *V, int size_1, int size_2);

/* For norm */
double matrixNorm(const double *a, int size_1, int size_2);
double discrepancy(const double *A, const double *x, const double *b, int size);

/* For multiplication */
void multiplyMatrices(
    const double *a, int a_m, int a_n,
    const double *b, int b_m, int b_n,
    double *c);

// void multiplyMatricesQuadratic(const double *a, const double *b, double *c, int m);

/* for sum (result will be in massive a) */
void addMatrices(double *a, const double *b, int size_1, int size_2);

/* for subtraction (result will be in massive a) */
void subtractMatrices(double *a, const double *b, int size_1, int size_2);

/* For printing */
void printMatrix(const double *A, int size_1, int size_2, int max_size);

/* 
  Solving Linear Equation System 
  1 - success, 0 - fault
*/
int solveLES(double *A, int matrix_size, int block_size, double *B, int *X_indeces, pthread_barrier_t *barrier, 
  int t_amount, int t_index, int formula_type, mainElem *main_elems);

/* Get (i;j) block from matrix*/
void copyBlockMatrix(const double *A, double *V, int i, int j, int matrix_size, int block_size, int v_size_1, int v_size_2);

/* Write (i;j) block in matrix */
void writeBlockMatrix(double *A, const double *V, int i, int j, int matrix_size, int block_size, int v_size_1, int v_size_2);

/* Change blocks */
void changeBlocks(double *A, int matrix_size, int block_size, double *V1, double *V2, int v_size_1, int v_size_2, 
  int change_pos_1_i, int change_pos_1_j, int change_pos_2_i, int change_pos_2_j);

/* Copy block from vector */
void copyBlockVector(const double *B, double *V, int i, int block_size, int v_size);

/* Write block in vector */
void writeBlockVector(double *B, const double *V, int i, int block_size, int v_size);

/* 
  Find main block in i row and swap it with current diagonal block in i-step 
  1 - success, 0 - fault
*/

double RecipMatrixNorm(double *String, int matrix_size, int block_size, int j, double matrix_A_norm, double *V1, double *V3);

/* Synchronyzation */
void reduce_sum(int p, mainElem *main_elems, int *swap_index, int *recip_found);

/* thread_func */
void *thread_func(void* arg);

/* work with time */
double get_time();
double get_full_time();

/* Copy and write blocks in string */
void writeBlockString(double *String, const double *V, int j, int matrix_size, int block_size, int v_size_1, int v_size_2);
void copyBlockString(const double *String, double *V, int j, int matrix_size, int block_size, int v_size_1, int v_size_2);

/* Swap quadratic blocks in string */
void swapBlockString(double *String, double *V1, double *V2, int index_1, int index_2, int matrix_size, int block_size);