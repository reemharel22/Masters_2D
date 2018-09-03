#include "problem.h"

inline double calculate_C(double R_k,double R_l) {
    return R_k * R_l;
}


/****
 * After each time step, we update the time parameters and check if we should stop
 * the calculation.
*/
int update_time(double*dt, double*time_passed, double*time_stop, int*cycle){
    *cycle += 1;
    *time_passed += *dt;

    // check if we should stop
    if (*time_passed == *time_stop) {
        return 1;
    }
    return 0;
}

double do_timestep(double dt, double **volume, double **R, double **Z, double ***A) {
    int i,j;
    
    return 0.0;
}

void build_matrix_A() {

}