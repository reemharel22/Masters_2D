#include "utils.h"
#include "problem.h"

void parse_data_file(char *f_name);

int main (int argc, char* argv[]) {
    double *E_current, *E_old;
    double *R, *Z;
    int K_max, L_max; // or N, M.
    double R_k, R_l, lambda, sigma;
    double Rho_1, Rho_2, Rho_3, Rho_4, C;
    E_current = malloc_1d( K_max * L_max );
    E_old     = malloc_1d( K_max * L_max );
    R         = malloc_1d( K_max * L_max );
    Z         = malloc_1d( K_max * L_max );


    // set up initial condition
    initial_conditions( K_max * L_max , E_old, R, Z);
    return 0;
}
