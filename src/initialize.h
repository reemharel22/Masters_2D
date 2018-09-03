#ifndef initialize_H_
#define initialize_H_

extern void init(int *K_max, int *L_max, int *KC_max,int *LC_max, 
double *dt, double *t0, int *cycle, double *time_passed,double *time_stop,
double ***R,double ***Z, double ***E_current, double ***E_old, double****A, double***volume);

void init_mesh_Kershaw1(int K, int L, double **R, double **Z);

void clean_prog(int K_max, int KC_max, int L_max, double **R, double **Z,
double **E_current, double **E_old, double***A,double**volume);

#endif