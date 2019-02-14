#include "problem.h"
#include "utils.h"
#include "math.h"
#include "calculations.h"
/**
 * @file problem.c
 */

/**
 * @brief End of time step, update time params.
 * @param t The Time struct.
 * checks if we got to the end time of the execution.
 * updates the cycle and the time passed.
 * aswell as calculates the new value of T
 */
int update_time(Time *t, Quantity *T) {
    double ** temp_curr = T->current;
    double ** temp_prev = T->prev;
    int X = T->KC_max;
    int Y = T->LC_max;
    double T1, T2, max_T = 0, min_T = 0, delta_temp, dt_tag, tmp = 0;
    int i, j;
    
    if (t->time_passed >= t->time_stop) {
        return 1;
    }

    for (i = 1; i < X; i++) {
        for (j = 1; j < Y; j++) {
            T1 = T->current[i][j];
            if (T1 > max_T) {
                max_T = T1;
            }
        }
    }
    
    min_T = (190 * 11605) * 1E-3;
    for (i = 1; i < X; i++) {
        for (j = 1; j < Y; j++) {
            T1 = T->current[i][j];
            T2 = T->prev[i][j];
            delta_temp = fabs(T2 - T1) / (T2 + min_T);
            if (delta_temp > tmp) {
                tmp = delta_temp;
            }
        }
    }
    dt_tag = t->dt_factor * t->dt / tmp;
    t->dt = minimum(dt_tag, 1.1*t->dt);
    if (t->dt > t->dt_max) {
        t->dt = t->dt_max;
    }

    t->cycle += 1;
    t->time_passed += t->dt;
    if (t->cycle % 1000 == 0) {
        printf("Ended time step number %d\n",t->cycle);
    }

    return 0;
}

/**
 * @brief the main time step.
 * 1. First we calculate the physical properties: Opacity, Diffusion Coefficient and the heat capacity.
 *
 *  
 */
void do_timestep(Problem *p) {
    int i, j;
    double ***A, **b;
    calculate_opacity(&p->opacity, &p->rho, &p->temp, &p->mats);
    calculate_diffusion_coefficient(&p->diff_coeff, &p->opacity, &p->constants, 0);
    calculate_heatcapacity(&p->heat_cap, &p->rho, &p->temp, &p->mats);

    apply_boundary(&p->energy.current,p->boundary_type, p->energy.KC_max, p->energy.LC_max);

    A = build_matrix_A(&p->coor, &p->vol, &p->diff_coeff, p->time.dt);
    b = build_b_vector(&p->energy, &p->temp, &p->opacity, &p->constants, p->time.dt);
    jacobi_method_naive(1000, p->energy.KC_max, p->energy.LC_max, 1e-10, A, p->energy.current, b);

    free_3d(A, p->energy.KC_max, p->energy.LC_max);
    free_2d(b, p->energy.LC_max);
    return;
}

/**
 * @brief Builds the matrix A by Kershaw scheme.
 * @param coor  the coordinates struct.
 * @param vol   the volume struct.
 * @param diff  the diffusion coefficient struct.
 * @param dt    the size of timestep.
 * 1. Calculates RK_KL RL_KL Z_KL Z_LK by equation 15c and 15d.
 * 2. Calculates the big SIGMA and big LAMBDA, although they are named sigma, lambda
 *    they are still big ones. (by equation 15a 15b).
 * 3. Calculates Rho[1~4] by equation 17e
 * 4. Calculates small sigma and small lambda.
 * 5. Calculates the matrix A by 17a-17f.
 */
double ***build_matrix_A(Coordinate *coor, Data *vol,Data *diff, double dt) {
    int i = 0,j = 0,k;

    int K_max = coor->K_max, L_max = coor->L_max, KC_max = vol->KC_max, LC_max = vol->LC_max;
        
    double **volume = vol->values;
    double **R = coor->R;
    double **Z = coor->Z;
    double **D = diff->values;

    double RK_KL[K_max][L_max];
    double RL_KL[K_max][L_max];
    double ZK_KL[K_max][L_max];
    double ZL_KL[K_max][L_max];
    
    double rho1[K_max - 1][L_max - 1];
    double rho2[K_max - 1][L_max - 1];
    double rho3[K_max - 1][L_max - 1];
    double rho4[K_max - 1][L_max - 1];

    double sigma[K_max][L_max], lambda[K_max][L_max];
    double jacobi[K_max][L_max];
    double ***A = malloc_3d(K_max, L_max, 10);
    
    //#pragma omp parallel for default(shared)
    for(i = 0; i < K_max - 1; i++){
        for (j = 0; j < L_max - 1; j++) {
            RK_KL[i][j] = (R[i + 1][j + 1] - R[i][j] + R[i + 1][j] - R[i][j + 1]) / 2.0;
            RL_KL[i][j] = (R[i + 1][j + 1] - R[i][j] - R[i + 1][j] + R[i][j + 1]) / 2.0;

            ZK_KL[i][j] = (Z[i + 1][j + 1] - Z[i][j] + Z[i + 1][j] - Z[i][j + 1]) / 2.0;
            ZL_KL[i][j] = (Z[i + 1][j + 1] - Z[i][j] - Z[i + 1][j] + Z[i][j + 1]) / 2.0;
            
            //calclulate jacobian..
            jacobi[i][j] = RK_KL[i][j] * ZL_KL[i][j] - RL_KL[i][j] * ZK_KL[i][j];
        }
    }
    
    // these sizes are cell quantities, therefore we go KC, LC
    // we do not go throught the last imaginary cells, ofc.
    // #pragma omp parallel for default(shared)
    for(i = 0; i < KC_max - 1; i++){
        for (j = 0; j < LC_max - 1; j++) {
            //the reason we take jacob(i,j) is because we already made sigma(i,j)
            // be the same index as jacob(i,j), only R is shifted. see (15)
            sigma[i][j] = (R[i + 1][j + 1] + R[i + 1][j]) 
            / ( D[i][j] * jacobi[i][j] + D[i + 1][j] * jacobi[i + 1][j]);

            lambda[i][j] = (R[i + 1][j + 1] + R[i][j + 1]) 
            / ( D[i][j] * jacobi[i][j] + D[i][j + 1] * jacobi[i][j + 1]);
            
        }
    }
    
    //maybe merge with the next loop for speed..
    //calculate the rhos
    double C_KL;
    double sqrt_sigma_ij,sqrt_sigma_i1j, sqrt_lambda_ij, sqrt_lambda_ij1;
    k = 0;
    for(i = 0; i < K_max - 1; i++){
        j = 0;
        sqrt_lambda_ij1 = sqrt_lambda_ij  = sqrt( lambda[i][j]);

        for (j = 0; j < L_max - 1; j++) {
            sqrt_sigma_ij   = sqrt( sigma[i][j]); 
            sqrt_sigma_i1j  = sqrt( sigma[k][j] );
            sqrt_lambda_ij  = sqrt( lambda[i][j]);
            C_KL            = RK_KL[i][j] + RL_KL[i][j];
            rho1[i][j]      = sqrt_sigma_ij  * sqrt_lambda_ij  * C_KL; 
            rho2[i][j]      = sqrt_sigma_i1j * sqrt_lambda_ij1 * C_KL; 
            rho3[i][j]      = sqrt_sigma_i1j * sqrt_lambda_ij  * C_KL;
            rho4[i][j]      = sqrt_sigma_ij  * sqrt_lambda_ij1 * C_KL;
            sqrt_lambda_ij1 = sqrt_lambda_ij;
        }
        k = i;
    }

    //calculate the small sigma and small lambda
    double r1,r2,r3,r4;
    double z1,z2,z3,z4;
    for(i = 0; i < KC_max - 1; i++){
        for (j = 0; j < LC_max - 1; j++) {
            r1 = RL_KL[i][j]     * RL_KL[i][j] ;
            r2 = RL_KL[i + 1][j] * RL_KL[i+1][j];
            r3 = RK_KL[i][j]     * RK_KL[i][j];
            r4 = RK_KL[i][j + 1] * RK_KL[i][j + 1];
            z1 = ZL_KL[i][j]     * ZL_KL[i][j];
            z2 = ZL_KL[i + 1][j] * ZL_KL[i + 1][j];
            z3 = ZK_KL[i][j]     * ZK_KL[i][j];
            z4 = ZK_KL[i][j + 1] * ZK_KL[i][j + 1];

            sigma[i][j] = sigma[i][j] * ( r1 + z1  + r2 + z2)  / 2.0;
            lambda[i][j] = lambda[i][j] * ( r3 + z3 + r4 + z4) / 2.0;
        }
    }

    // the matrix is represnted in the following manner:
    // 6 7 8
    // 3 4 5
    // 0 1 2    9 is the solution.. i guess
    double vol_1;
    double a5,a6,a7,a8;
    int p = 0;
    k = 0;
    for(i = 0; i < K_max - 1; i++){
        //k is equivilent to i - 1. instead of MAX(i,0)...
        for (j = 0; j < L_max - 1; j++) {
            vol_1 = 1.0 / volume[i][j];
            A[i][j][4] = 1.0 - dt*vol_1*( -sigma[i][j] - sigma[k][j] - lambda[i][j] - lambda[i][j - 1]
                         + 0.5* (rho1[i][j] + rho2[i][j] - rho3[i][j] - rho4[i][j]) );
            
            a8 = -0.25* dt * (rho3[i + 1][j] + rho4[i][j + 1]);
            a7 = dt * (lambda[i][j] - 0.25*(rho1[i][j] + rho2[i][j + 1] - rho3[i][j] - rho4[i][j + 1]));  
            a6 = 0.25 *dt * (rho1[k][j] + rho2[i][j + 1]);
            a5 = dt *(sigma[i][j] - 0.25 * (rho1[i][j] + rho2[i][j + 1] - rho3[i][j] - rho4[i][j + 1]));
            
            A[i][j][8]         = vol_1 * a8; // equation (17d)
            A[i][j][7]         = vol_1 * a7; // equation (17c)
            A[i][j][6]         = vol_1 * a6; // equation (17e)
            A[i][j][5]         = vol_1 * a5; // equation (17b)
           
            A[k][p][0]         = (1.0/volume[k][p])     * a8;
            A[i][p][1]         = (1.0/volume[i][p])     * a7;
            A[i + 1][p][2]     = (1.0/volume[i + 1][p]) * a6;
            A[k][j][3]         = (1.0/volume[k][j])     * a5;
            p = j;
        }
        p = 0;
        k = i;
    }

    return A;
}

/**
 * @brief The first boundary condition we apply is the "1d" effective one.
 * In which the bottom and top walls are "reflective", the left corner gets 190ev and right corner is 0.
 * 
 * 
*/
void apply_boundary(double **data, int boundary_type, int n, int m) {
    int i, j;
    //TODO check what does it mean n*n*RK?????
    //TODO tell shay, as oppose to the paper. i calculate it like "them"
    //i.e boundary condition is abit different


    // TOP WALL - reflective
    // i and j = m - 1
    for (i = 0; i < n; i++) {
        data[i][m - 1]     = data[i][m - 2];
    }

    // Bottom WALL - reflective
    // i and j = 0
    for (i = 0; i < n; i++) {
        data[i][0]     = data[i][1];
    }

    //RIGHT WALL - leapes
    for (j = 0; j < m; j++) {
        data[n - 2][j] = 0;
    }
    
    //LEFT WALL - 190 ev
    for (j = 0; j < m; j++) {
        data[0][j] = 190 * 11605;
    }
}

/**
 * @brief Prepares the solution vector b, Ax=b
 * energy - is the prev energy i.e E(r, t = n), current = n + 1 
 * 
*/
double **build_b_vector(Quantity *E, Quantity *T, Data *opac, Constants *consts, double dt) {
    int i, j;
    int X = T->KC_max;
    int Y = T->LC_max;
    double c = consts->c_light;
    double **b = malloc_2d(X, Y);
    double **energy = E->prev;
    double **temp = T->prev;
    double **opacity =  opac->values;
    for (i = 0; i < X; i++) {
        for (j = 0; j < Y; j++) {
            b[i][j] = energy[i][j] + opacity[i][j] * dt * c * temp[i][j];
        }
    }
    return b;
}

