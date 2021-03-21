#include <sys/time.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "utils.h"

extern const int VERBOSE = 1;
extern const double HEV = 1160500;
extern const double PI = 3.14159265359;

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

void check_nan_2d(double **arr, int nx, int ny, char* str){
    int i, j; 
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            if (arr[i][j] != arr[i][j]) {
                printf("Found NaN in %s!.\t  Position: with i: %d and j: %d\n",str,  i, j);
                exit(1);
            }
        }
    }
    // printf("NaN not found!\n");
}

void check_nan_3d(double ***arr, int nx, int ny, int nz) {
    int i, j,k; 
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nz; k++) {
                if (arr[i][j][k] != arr[i][j][k]) {
                    printf("Found NaN 3d with i: %d and j: %d\n", i,j);
                    exit(1);
                    return;
                }
            }
        }
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

void print_matrix_row(int nx, int ny, double**val) {
    int i,j;
    for( i = 0; i < nx; i++) {
        for( j = 0; j < ny; j++){
            printf("%lf\n",val[i][j]);
        }
        printf("\nend row: %d\n",i);
    }
} 

void print_matrix_col(int nx, int ny, double**val) {
    int i,j;
    
    for( j = 0; j < ny; j++){
        for( i = 0; i < nx; i++) {
            printf("%lf\n",val[i][j]);
        }
        printf("end col: %d\n",j);
    }
    
} 

void print(double arr) {
    printf("%16e\n", arr);

}

void print3(double **arr, int i, int j){
    printf("%10e\t %d \t %d\n",arr[i][j],i,j);
}


void print_2d(double **arr, int nx, int ny) {
    int i,j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            printf("%10e\t", arr[i][j]);
        }
        printf("\n");
    }
}

void print_3d(double ***arr, int nx, int ny, int Q) {
    int i, j, k;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < Q; k++) {
                printf("%10e\t", arr[i][j][k]);
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
    return 0;
}

void mesh_square_volume(double **volume, double **X, double **Y, 
int n, int m) {
    int i,j;
    double x1,x2,y1,y2;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            y2 = Y[i][j+1];
            y1 = Y[i][j];
            x2 = X[i+1][j];
            x1 = X[i][j];
            // if (n == 3){ // 1d problem on y axis
            //     volume[i][j] = y2 - y1;
            // }
            // else if (m == 3) { // 1d problem on x axis
            //     volume[i][j] = x2 - x1;
            // }
            // else {
                volume[i][j] = fabs((y2 - y1) * PI*(x2*x2 - x1*x1)) / (2 * PI);
                // volume[i][j] =  (y2-y1) * (x2 - x1) ;

            // }
        }
    }

}

void solveTriagonal(int N, double * solve, double *L,double *U, double* mainD) {
    int i;
    U[0] = U[0] / mainD[0];
    solve[0] = solve[0] / mainD[0];

    /* loop from 1 to X - 1 inclusive, performing the forward sweep */
    for (i = 1; i < N; i++) {
        const double m = 1.0 / (mainD[i] - L[i] * U[i - 1]);
        U[i] = U[i] * m;
        solve[i] = (solve[i] - (L[i] * (solve[i - 1]))) * m;
    }

    /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
    for (i = N - 2; i >=0 ; i--) {
        solve[i] -= U[i] * solve[i + 1];
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
void gauss_seidel_method(int max_iter, int nx, int ny, double epsilon, double ***A, double **x, double **b) {
    int i, j;
    double sum;
    int iter;
    double aii_r;
    for (iter = 0; iter < max_iter; iter++) {
        // n^2 is i to j...
        for (i = 0; i < nx; i++) {
            for (j = 0; j < nx; j++) {
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
void jacobi_method_naive(int max_iter, int nx, int ny, double epsilon, double ***A, double **x, double **b) {
    int i, j, iter;
    double sum;
    double aii_r;
    double ** x_prev;
    int im,jm;
    x_prev = malloc_2d(nx, ny);
    double **tmp;
    for (i = 0; i < nx; i++){
        for (j = 0; j < ny; j++) {
            x_prev[i][j] = x[i][j];
        }
    }
    // print_2d(b,nx,ny);
    // 6 7 8
    // 3 4 5
    // 0 1 2    9 is the solution.. i guess
    for (iter = 0; iter < max_iter; iter++) {
        for (i = 1; i < nx - 1; i++) {
            for (j = 1; j < ny - 1; j++) {
                sum = 0.0;
                im = maximum(i - 1, 0);
                jm = maximum(j - 1, 0);
                //let's write it explicitly Maybe should be A[i][j]
                sum += A[i][j][0] * x_prev[i + 1][jm]; // A[i][j][0] is left bottom.
                sum += A[i][j][1] * x_prev[i + 1][j];     // A[i][j][1] is bottom.
                sum += A[i][j][2] * x_prev[i + 1][j + 1]; // A[i][j][2] is right bottom.
                sum += A[i][j][3] * x_prev[i][jm];     // A[i][j][3] is left.
                sum += A[i][j][5] * x_prev[i][j + 1];     // A[i][j][5] is right.
                sum += A[i][j][6] * x_prev[im][jm]; // A[i][j][6] is left top.
                sum += A[i][j][7] * x_prev[im][j];     // A[i][j][7] is top.
                sum += A[i][j][8] * x_prev[im][j + 1]; // A[i][j][8] is top right.
                aii_r = 1.0 / A[i][j][4];
                x[i][j] = aii_r * (b[i][j] - sum);
            }
        }
        if (converge(nx, ny, epsilon, x_prev, x)) {
            // print_2d(x_prev,nx,ny);
            // printf("\n\n");
            // print_2d(x,nx, ny);
            #pragma omp parallel for collapse(2)
            for (i = 1; i < nx - 1; i++){
                for (j = 1; j < ny - 1; j++) {
                    x_prev[i][j] = x[i][j];
                }
            }
             free_2d(x_prev, nx);
             printf("Converged after %d\n", iter);
             return;
         }
        #pragma omp parallel for collapse(2)
            for (i = 1; i < nx-1; i++){
            for (j = 1; j < ny-1; j++) {
                    x_prev[i][j] = x[i][j];
                }
            }
        }
    free_2d(x_prev, nx );
}

void check_boundary_condition(double x, double bc_value) {
    int i, j;
    
    if (abs(x - bc_value) / abs(x + bc_value) > 0.01) 
    {
        printf("Boundary condition doesn't apply. \nValue: %10e whereas BC is :%10e\n", x, bc_value);
        // exit(1);
    }
    printf("Boundary condition works!\n");
}

void mat_mul(double ***A, double ** x, double **b, double epsilon, int nx, int ny) {
    int i, j, k = 0;
    double sum;
    // solves the Ax=b equation where x is the solution from the jacobi method... obviously :)
    double **vec = malloc_2d(nx, ny);
    for (i = 1 ; i < nx - 1; i++ ) {
        for (j = 1; j < ny - 1; j ++) {
            sum = 0.0;
            sum += A[i][j][0] * x[i + 1][j - 1]; // A[i][j][0] is left bottom.
            sum += A[i][j][1] * x[i + 1][j];     // A[i][j][1] is bottom.
            sum += A[i][j][2] * x[i + 1][j + 1]; // A[i][j][2] is right bottom.
            sum += A[i][j][3] * x[i][j - 1];     // A[i][j][3] is left.
            sum += A[i][j][5] * x[i][j + 1];     // A[i][j][5] is right.
            sum += A[i][j][6] * x[i - 1][j - 1]; // A[i][j][6] is left top.
            sum += A[i][j][7] * x[i - 1][j];     // A[i][j][7] is top.
            sum += A[i][j][8] * x[i - 1][j + 1]; // A[i][j][8] is top right.
            sum += A[i][j][4] * x[i][j];
            vec[i][j] = sum;
        }
    }
    if (converge(nx,ny,epsilon, vec, b)) {
        printf("Ax=b worked!\n");
    }
    else{
        print_2d(vec, nx, ny);
        printf("\n\n");
        print_2d(b, nx, ny);
        printf("didnt solve correctly\n");
        exit(1);
    }
}


void read_file(char * f_name, double * x, int nx) {
    int i;
    FILE *fp;
    fp = fopen(f_name, "r");
    fscanf(fp, "%lf", &x[0]);
    
    for (size_t i = 0; i < nx - 1; i++)
    {
        fscanf(fp, "%lf", &x[i]);
        // print(x[i]);
    }
    
}

void check_monotoic_up(double **x, int nx, int ny) {
    int i,j;
    for (i = 2; i < nx - 1; i++){
        for (j = 2; j < ny - 1; j++) {
            if ((1.001*x[i][j] ) < x[i][j +1]) {
                printf("Not monotic ole %d %d\n", i, j);
                print(x[i][j]);
                print(x[i+1][j]);
                exit(1);
            }
        }
    }
}

int converge(int nx, int ny, double epsilon, double **x_prev, double** x_current) {
    int i, j;
    double diff, sqdiff = 0.0;
    for (i = 1; i < nx - 1; i++) {
        for (j = 1; j < ny - 1; j++) {
            diff = (1 - x_current[i][j] / x_prev[i][j]);
            // sqdiff = diff * diff;
            // printf ("diff ");
            // print(fabs(diff));
            // print(epsilon);
            // printf("%10e %10e %10e %10e \n", x_prev[i][j], x_current[i][j], diff, epsilon);
            if (fabs(diff) > epsilon) {
                return 0;
            }
        }
    }
    return 1;
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