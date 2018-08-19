#include <sys/time.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double get_time() {
   struct timeval stop, start;
    gettimeofday(&start, 0);
    return start.tv_usec;
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

void write_to_file(char* f_name, double * values, int n) {
    int i;
    FILE *f = fopen(f_name, "a+");
    
    for(i = 0; i < n; i++)
    {
        fprintf(f, "%lf ", values[i]);
    }
    
}