#ifndef problem_H_
#define problem_H_
#include "structs.h"
void initial_conditions(int, double*, double*, double *);
void apply_boundary_condition(double*);
inline double calculate_C(double,double);
void do_timestep(Problem *p);
int update_time(Time *t);
double ***build_matrix_A(Coordinate *coor, Volume *vol,Diff_Coeff *diff, double dt);
#endif
