#include <sys/time.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "utils.h"
#include <omp.h>

inline double get_time() {
  return omp_get_wtime();
}

inline double avg(double x, double y) {
    return (x + y) * 0.5;
}

inline double avg_harmonic(double x, double y) {
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
    double * ptr =  (double*) malloc(sizeof(double) * n);
    return ptr;
}

double** malloc_2d(int n, int m) {
    int i = 0;
    double ** arr = (double **) malloc(n * sizeof(double *));
    for( i = 0; i < n; i++) {
        arr[i] = malloc_1d(m);
    }
    return arr;
}

double*** malloc_3d(int n, int m, int q) {
    int i = 0;
    double *** arr = (double ***) malloc(n * sizeof(double **));
    for( i = 0; i < n; i++) {
        arr[i] = malloc_2d(m,q);
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
            fprintf(f, "%lf %lf\n",X[i][j],Y[i][j]);
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

inline double square_volume(double x1, double x2, double y1, double y2) {
    return (y2-y1) * (x2-x1);
}

void mesh_square_volume(double **volume, double **X, double **Y, 
int n, int m) {
    int i,j;
    for ( i = 0; i < n-1; i++) {
        for (j = 0; j < m-1; j++) {
            volume[i][j] = square_volume(X[i][j],X[i][j + 1],Y[i][j],Y[i + 1][j]);           
        }
    }
}
