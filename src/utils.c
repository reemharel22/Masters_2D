#include <sys/time.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "utils.h"

double get_time() {
  return omp_get_wtime();
}

double avg(double x, double y) {
    return (x + y) * 0.5;
}

double avg_harmonic(double x, double y) {
    if (x == 0) {
        printf("warning utils.c - x is zero, replacing it with 1e-15\n");
        x = 1-15;
    }
    if (y == 0) {
        printf("warning utils.c - y is zero, replacing it with 1e-15\n");
        y = 1-15;
    }
   return (2.0 * x * y)/(x + y);
}

double* malloc_1d(int n) {
    double * ptr =  malloc(sizeof(double) * n);
    return ptr;
}

double** malloc_2d(int n, int m) {
    int i = 0;
    int j = 0;
    double ** arr = (double **) malloc(n * sizeof(double *));
    for( i = 0; i < n; i++) {
        arr[i] = malloc_1d(m);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            arr[i][j] = 0.0;
        }
    }
    return arr;
}

double*** malloc_3d(int n, int m, int q) {
    int i = 0;
    int j = 0;
    int k = 0;
    double *** arr = (double ***) malloc(n * sizeof(double **));
    for( i = 0; i < n; i++) {
        arr[i] = malloc_2d(m,q);
    }
     for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            for (k = 0; k < q; k++) {
                arr[i][j][k] = 0.0;
            }
        }
    }
    return arr;
}

void write_to_file(char* f_name, double * values, int n) {
    int i;
    FILE *f = fopen(f_name, "a+");
    
    for(i = 0; i < n; i++)
    {
        fprintf(f, "%lf ", values[i]);
    }
    
}

void write_to_file_mesh(char * f_name, double **X, double **
Y, int n, int m) {
    int i,j;
    FILE *f = fopen(f_name, "w+");
    for (i = 0 ; i < n; i++ ) {
        for (j = 0; j < n; j++) {
            fprintf(f, "%lf %lf\n", X[i][j], Y[i][j]);
        }
    }
    fclose(f);
}

void print_matrix_row(int K, int L, double**val) {
    int i,j;
    for( i = 0; i < K; i++) {
        for( j = 0; j < L; j++){
            printf("%lf\n",val[i][j]);
        }
        printf("\nend row: %d\n",i);
    }
} 

void print_matrix_col(int K, int L, double**val) {
    int i,j;
    
    for( j = 0; j < L; j++){
        for( i = 0; i < K; i++) {
            printf("%lf\n",val[i][j]);
        }
        printf("end col: %d\n",j);
    }
    
} 


void print_2d(double **arr, int K, int L) {
    int i,j;
    for (i = 0; i < K; i++) {
        for (j = 0; j < L; j++) {
            printf("%10e\t", arr[i][j]);
        }
        printf("\n");
    }
}

void print_3d(int K, int L, int Q, double ***arr) {
    int i, j, k;
    for (i = 0; i < K; i++) {
        for (j = 0; j < L; j++) {
            for (k = 0; k < Q; k++) {
                printf("%10e\t", arr[i][j]);
            }
            
        }
        printf("\n");
    }
}


void free_1d(double *ptr) {
    free(ptr);
}

void free_2d(double **ptr, int n) {
    int i = 0;
    for ( i = 0; i < n; i ++) {
        free(ptr[i]);
    }
    free(ptr);
}

void free_3d(double ***ptr, int n, int m) {
    int i = 0,j;
    for ( i = 0; i < n; i ++) {
        for (j = 0; j < m; j++)
            free(ptr[i][j]);
        free(ptr[i]);
    }
    free(ptr);
}

double square_volume(double x1, double x2, double y1, double y2) {
    return (y2-y1) * (x2-x1);
}

void mesh_square_volume(double **volume, double **X, double **Y, 
int n, int m) {
    int i,j;
    for ( i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            volume[i][j] = square_volume(X[i][j],X[i][j + 1],Y[i][j],Y[i + 1][j]);         
        }
    }
}

/**
 * @brief Iterative matrix Ax=b solver in Gauss-Seidel method.
 * iter = maximum number of iterations
 * n = size of matrix, and vectors
 * epsilon accuaracy requiremnt
 * A the matrix
 * x pointer to the the solved vector.
 * b pointer to the vector
 * input is a pointer to the matrix
 */
void gauss_seidel_method(int max_iter, int K, int L, double epsilon, double ***A, double **x, double **b) {
    int i, j;
    double sum;
    int iter;
    double aii_r;
    for (iter = 0; iter < max_iter; iter++) {
        // n^2 is i to j...
        for (i = 0; i < K; i++) {
            for (j = 0; j < K; j++) {
                sum = 0.0;
                //for (k = 0; k < 10; k++) {
                //let's write it explicitly
                sum += A[i][j][0] * x[i + 1][j - 1]; // A[i][j][0] is left bottom.
                sum += A[i][j][1] * x[i + 1][j];     // A[i][j][1] is bottom.
                sum += A[i][j][2] * x[i + 1][j + 1]; // A[i][j][2] is right bottom.
                sum += A[i][j][3] * x[i][j - 1];     // A[i][j][3] is left.
                sum += A[i][j][5] * x[i][j + 1];     // A[i][j][5] is right.
                sum += A[i][j][6] * x[i - 1][j - 1]; // A[i][j][6] is left top.
                sum += A[i][j][7] * x[i - 1][j];     // A[i][j][7] is top.
                sum += A[i][j][8] * x[i - 1][j + 1]; // A[i][j][8] is top right.
                //}
                aii_r = 1.0 / A[i][j][4];
                x[i][j] = aii_r * (b[i][j] - sum);
            }
        }
    }
}


/**
 * @brief Iterative matrix Ax=b solver in Gauss-Seidel method suitable for only kershaw scheme.
 * iter = maximum number of iterations
 * k * l = size of matrix, and vectors
 * epsilon accuaracy requiremnt
 * A the matrix size 10 * k * l
 * x pointer to the the solution vector.
 * b pointer to the vector
 * input is a pointer to the matrix
 * 
 */
void jacobi_method_naive(int max_iter, int K, int L, double epsilon, double ***A, double **x, double **b) {
    int i, j, iter;
    double sum;
    double aii_r;
    double ** x_prev;
    x_prev = malloc_2d(K + 1, L + 1);
    double **tmp;
    for (i = 0; i < K; i++){
        for (j = 0; j < L; j++) {
            x_prev[i][j] = x[i][j];
        }
    }
  //  print_2d(K, L, x_prev);
    for (iter = 0; iter < max_iter; iter++) {
        // n^2 is i to j...
        for (i = 1; i < K; i++) {
            sum = 0.0;
            for (j = 1; j < L; j++) {
                //for (k = 0; k < 10; k++) {
                //let's write it explicitly Maybe should be A[i][j]
                sum += A[i][j][0] * x_prev[i + 1][j - 1]; // A[i][j][0] is left bottom.
                sum += A[i][j][1] * x_prev[i + 1][j];     // A[i][j][1] is bottom.
                sum += A[i][j][2] * x_prev[i + 1][j + 1]; // A[i][j][2] is right bottom.
                sum += A[i][j][3] * x_prev[i][j - 1];     // A[i][j][3] is left.
                sum += A[i][j][5] * x_prev[i][j + 1];     // A[i][j][5] is right.
                sum += A[i][j][6] * x_prev[i - 1][j - 1]; // A[i][j][6] is left top.
                sum += A[i][j][7] * x_prev[i - 1][j];     // A[i][j][7] is top.
                sum += A[i][j][8] * x_prev[i - 1][j + 1]; // A[i][j][8] is top right.
                aii_r = 1.0 / A[i - 1][j - 1][4];
    //            printf("%10e\n", A[i - 1][j - 1][1]);
                x[i][j] = aii_r * (b[i - 1][j - 1] - sum);
            }
        }
        if (converge(K, L, epsilon, x_prev, x)) {
            free_2d(x_prev, K);
            return;
        }
        tmp = x_prev;
        x_prev = x;
        x = tmp;
    }
    free_2d(x_prev, K);
}

int converge(int K, int L, double epsilon, double **x_prev, double** x_current) {
    int i, j;
    double diff,sqdiff = 0.0;
    for (i = 0; i < K; i++) {
        for (j = 0; j < L; j++) {
            diff = x_prev[i][j] - x_current[i][j];
            sqdiff += diff * diff;
        }
    }
    if (sqdiff > epsilon) {
        return 0;
    } else {
        return 1;
    }
}

double minimum(double a, double b) {
    if (a < b) {
        return a;
    }
    return b;
}

double maximum(double a, double b) {
    if (a > b) {
        return a;
    }
    return b;
}

int int_reader(char *str, int len) {
    char delim[] = " ";
    char *ptr = strtok(str, delim);
    ptr = strtok(NULL, delim);
    int dble = atoi(ptr);
    return dble;
}

int* int_array_reader(char *str, int len) {
    char delim[] = " ";
    len = -1;
    char * tmp[50];
    strcpy(tmp, str);
    char *ptr = strtok(tmp, delim);
    while (ptr != NULL) {
        len ++;
        ptr = strtok (NULL, " ");
    }
    int * arr = (int*) malloc(sizeof(int) * len);
    char *ptr1 = strtok(str, delim);
    int i = 0;
    for ( i = 0; i < len; i++) {
        ptr1 = strtok (NULL, " ");
        arr[i] = atoi(ptr1);
    }
    return arr;
}

double double_reader(char *str, int len) {
    char delim[] = " ";
    char *ptr = strtok(str, delim);
    ptr = strtok(NULL, delim);
    double dble = atof(ptr);
    return dble;
}

double* double_array_reader(char *str, int len) {
    char delim[] = " ";
    len = -1;
    char * tmp[50];
    strcpy(tmp, str);
    char *ptr = strtok(tmp, delim);
    while (ptr != NULL) {
        len ++;
        ptr = strtok (NULL, " ");
    }
    double * arr = (double*) malloc(sizeof(double) * len);
    char *ptr1 = strtok(str, delim);
    int i = 0;
    for ( i = 0; i < len; i++) {
        ptr1 = strtok (NULL, " ");
        arr[i] = atof(ptr1);
    }
    return arr;
}