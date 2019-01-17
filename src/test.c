#include <sys/time.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "utils.h"

void prep_b(double **b, int K, int L) {
    int i,j;
    int t1;
    int t2;
    for (i = 0; i < K; i++) {
        for (j = 0; j < L; j++) {
            t1 = i + j;
            t2 = L + K;
            b[i][j] = t2 / (t1 + 1);
        }
    }
}


void prep_A(double ***A, int K, int L) {
    int i, j;
    for (i = 0; i < K; i++) {
        for (j = 0; j < L; j++) {
            A[i][j][4] = 1.0;
        }
    }
    A[0][0][1] = -1;
    A[0][0][2] = 1;
    A[0][0][4] = 1;
    A[0][0][5] = 1.2;

    //A[0][1][0] = 0.5;
    A[0][1][1] = 0.5;
    A[0][1][2] = -0.5;
    A[0][1][3] = -1;
    A[0][1][5] = -2;
    A[0][1][6] = 1;
    A[0][1][7] = 1.2;
    
    A[0][2][1] = 5;
    A[0][2][2] = 1;
    A[0][2][3] = 0.5;
    A[0][2][5] = 0.5;
    A[0][2][6] = 1;
    A[0][2][7] = -2;

    
    A[0][3][2] = 1;
    A[0][3][3] = 5;
    A[0][3][5] = 1.7;
    A[0][3][6] = 1;
    A[0][2][7] = 0.5;

    A[1][0][2] = 1;
    A[1][0][5] = 2;
    A[1][0][6] = 1;
    A[1][0][7] = 1.7;
    
    A[1][1][1] = 2.2;
    A[1][1][2] = 1;
    A[1][1][6] = 1;
    A[1][1][7] = 2;
    
    A[1][2][1] = 1.2;
    A[1][2][2] = 1;
    A[1][2][3] = 2.2;
    A[1][2][5] = 0;
    A[1][2][6] = 1;
    A[1][2][7] = 0;

    A[1][3][1] = 0;
    A[1][3][2] = 1;
    A[1][3][3] = 1.2;
    A[1][3][5] = 0;
    A[1][3][6] = 1;
    A[1][3][7] = 0;

    A[2][0][1] = 1.5;
    A[2][0][2] = 1;
    A[2][0][5] = 0.5;
    A[2][0][6] = 1;
    
    A[2][1][2] = 1;
    A[2][1][3] = 1.5;
    A[2][1][6] = 1;
    A[2][1][7] = 0.5;

    A[2][2][1] = 2;
    A[2][2][2] = 1;
    A[2][2][5] = 0;
    A[2][2][6] = 1;

    A[2][3][3] = 2;
    A[2][3][6] = 1;
}

int main() {
    int K = 3;
    int L = 4;
    int max_iter = 1000;
    double ***A;
    double **x;
    double **b;
    double epsilon = 1E-4;
    x = malloc_2d(K + 2, L + 2);
    A = malloc_3d(K, L, 10);
    b = malloc_2d(K + 2, L + 2);
    //p/rep_b(b, K, L);
    //prep_A(A, K, L);
    //jacobi_method_naive(max_iter, K + 1, L + 1, epsilon, A, x, b);
    //print_2d(K + 2, L + 2, x); 
    int i = 0;
    #pragma omp parallel
    {
        #pragma omp for
        for (i = 0; i < K * L; i++) {
            printf("%d\n", omp_get_thread_num());
        }
    }
    return 0;
}