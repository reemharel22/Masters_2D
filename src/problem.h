#ifndef problem_H_
#define problem_H_
#include "structs.h"

void do_timestep(Problem *p);
int update_time(Time *t, Quantity *T);
void apply_boundary(double **data,int, int n, int m);
/**
 * @brief Builds the matrix A by Kershaw scheme.
 * @var Coordinate::coor  the coordinates struct.
 * @struct vol   the volume struct.
 * @struct diff  the diffusion coefficient struct.
 * @struct dt    the size of timestep.
 *           1. Calculates RK_KL RL_KL Z_KL Z_LK by equation 15c and 15d.
 *           2. Calculates the big SIGMA and big LAMBDA, although they are named sigma, lambda
 *              they are still big ones. (by equation 15a 15b).
 *           3. Calculates Rho[1~4] by equation 17e
 *           4. Calculates small sigma and small lambda.
 *           5. Calculates the matrix A by 17a-17f.
 */
double ***build_matrix_A(Coordinate *coor, Data *vol,Data *diff, double dt);
double **build_b_vector(Quantity *E, Quantity *T, Data *opac, Constants *consts, double dt);
void update_prev_values(Problem * p);
#endif
