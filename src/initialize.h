#ifndef initialize_H_
#define initialize_H_
#include "structs.h"
/*extern void init(int *K_max, int *L_max, int *KC_max,int *LC_max, 
double *dt, double *t0, int *cycle, double *time_passed,double *time_stop,
double ***R,double ***Z, double ***E_current, double ***E_old, double****A, double***volume);*/

void init(Problem*p, char*argv);
void init_mesh_Kershaw1(int K, int L, double **R, double **Z, double dx);
void init_density(Materials *mats, Data *density);
void clean_prog(Problem *p);
void init_datafile(Problem *p, char* f_name);
void init_materials_datafile(Material *m, int mat_number);
#endif
