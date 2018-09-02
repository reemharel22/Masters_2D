#ifndef utils_H_
#define utils_H_

double get_time();
inline double avg(double x, double y);
inline double avg_harmonic(double, double);
double* malloc_1d(int n);
double** malloc_2d(int n, int m);
void write_to_file(char* f_name, double * values, int n);

#endif