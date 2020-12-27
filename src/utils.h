#ifndef utils_H_
#define utils_H_

double get_time();
double avg(double x, double y);
double avg_harmonic(double, double);
double* malloc_1d(int n);
double** malloc_2d(int n, int m);
double*** malloc_3d(int n, int m,int);
void write_to_file(char* f_name, double * values, int n);
void write_to_file_mesh(char * f_name, double **X, double **Y, int n, int m);
void print_matrix_row(int K, int L, double**val) ;
void print_matrix_col(int K, int L, double**val) ;
void free_2d(double **ptr, int n);
void free_3d(double ***ptr, int n,int m);
void print_2d(double **arr, int K, int L);
void print_3d(double ***arr, int K, int L, int Q);
void free_1d(double *ptr);
double minimum(double, double);
double maximum(double, double);
double square_volume(double x1, double x2, double y1, double y2);
void mesh_square_volume(double **volume, double **X, double **Y, int n, int m);
int converge(int K, int L, double epsilon, double **x_prev, double** x_current);
void jacobi_method_naive(int max_iter, int K, int L , double epsilon, double ***A, double **x, double **b);
void gauss_seidel_method(int max_iter, int K, int L, double epsilon, double ***A, double **x, double **b);
int int_reader(char *str, int len);
int* int_array_reader(char *str, int len);
double double_reader(char *str, int len);
double* double_array_reader(char *str, int len);
void check_nan_2d(double **arr, int K, int L, char* s);
void check_nan_3d(double ***arr, int K, int L, int Q);
#endif