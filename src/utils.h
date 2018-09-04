#ifndef utils_H_
#define utils_H_

inline double get_time();
inline double avg(double x, double y);
inline double avg_harmonic(double, double);
double* malloc_1d(int n);
double** malloc_2d(int n, int m);
double*** malloc_3d(int n, int m,int);
void write_to_file(char* f_name, double * values, int n);
void write_to_file_mesh(char * f_name, double **X, double **Y, int n, int m);
void print_matrix_row(int K, int L, double**val) ;
void print_matrix_col(int K, int L, double**val) ;
void free_2d(double **ptr, int n);
void free_3d(double ***ptr, int n,int m);
void free_1d(double *ptr);
inline double square_volume(double x1, double x2, double y1, double y2);
void mesh_square_volume(double **volume, double **X, double **Y, int n, int m);
#endif