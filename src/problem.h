#ifndef problem_H_
#define problem_H_

void initial_conditions(int, double*, double*, double *);
void apply_boundary_condition(double*);
inline double calculate_C(double,double);
double do_timestep();
int update_time(double*dt, double*time_passed, double*time_stop, int*cycle);
void build_matrix_A();
#endif
