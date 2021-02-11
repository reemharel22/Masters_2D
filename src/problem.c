#include "problem.h"
#include "utils.h"
#include <stdio.h>
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
    int X = T->nx;
    int Y = T->ny;
    double T1, T2, max_T = 0, min_T = 0, delta_temp, dt_tag, tmp = 0;
    int i, j;
    
    if (t->time_passed >= t->time_stop) {
        return 1;
    }

    for (i = 1; i < X - 1; i++) {
        for (j = 1; j < Y - 1; j++) {
            T1 = T->current[i][j];
            if (T1 > max_T) {
                max_T = T1;
            }
        }
    }
    
    min_T = max_T * 1E-3;
    for (i = 1; i < X - 1; i++) {
        for (j = 1; j < Y - 1; j++) {
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

    return 0;
}

/**
 * @brief the main time step.
 * 1. First we calculate the physical properties: Opacity, Diffusion Coefficient and the heat capacity.
 * THE TEMPERATURE ISNT REALLY T BUT V...
 *  
 */
void do_timestep(Problem *p) {
    double ***A, **b;
    int i,j;
    int nx,ny;

    nx = p->temp->nx;
    ny = p->temp->ny;
    update_prev_values(p);

    //BC FOR OLSON
    for (j = 0; j < ny; j++){
        p->energy->prev[0][j] = pow(p->constants->TH,4) * p->constants->a_rad;
        p->temp->prev[0][j] = pow(p->constants->TH,4) * p->constants->a_rad;
    }

    if (VERBOSE)
        printf("Calculating Opacity\n");
    calculate_opacity(p->opacity, p->rho, p->temp, p->mats); //with prev

    if (VERBOSE)
        printf("Calculating The diffusion coefficient\n");
    calculate_diffusion_coefficient(p->diff_coeff, p->opacity, p->constants, 0);

    if (VERBOSE)
        printf("Calculating heat capacity\n");
    calculate_heat_capacity(p->heat_cap, p->rho, p->temp, p->mats);

    //apply boundary on the energy and temperature
    if (VERBOSE)
        printf("Building A and b matrix\n");
    //Calculating the energy.

    A = build_matrix_A(p->coor, p->vol, p->diff_coeff, p->time->dt, p->constants, p->energy);
    // print_3d(A, p->temp->nx, p->temp->ny, 10);
    b = build_b_vector(p->energy, p->temp, p->opacity,p->coor->R, p->constants, p->time->dt, p->time->cycle);
    // print_2d(b,p->energy->nx ,p->energy->ny );

    if (VERBOSE)
        printf("Solving the matrix\n");
    jacobi_method_naive(1000, p->energy->nx, p->energy->ny, 0.0001, A, p->energy->current, b);

    mat_mul(A, p->energy->current, b, 0.001, p->energy->nx, p->energy->ny);

    free_3d(A,p->energy->nx + 1, p->energy->ny + 1);
    check_nan_2d(p->energy->current, p->energy->nx ,p->energy->ny, "energy");
    if (VERBOSE)
        printf("Calculating Temperature\n");
    calculate_temperature(p->temp, p->energy, p->constants, p->opacity, p->heat_cap, p->time->dt);
        for (j = 0; j < ny-1; j++){
    // print_2d(p->energy->current, p->energy->nx ,p->energy->ny);

        // check_boundary_condition(p->energy->current[1][j], pow(p->constants->TH, 4) * p->constants->a_rad);
    }
    if (VERBOSE)
        printf("Updating the current prev values\n");
    check_nan_2d(p->temp->current, p->energy->nx ,p->energy->ny, "temperature");
    
    free_2d(b, p->energy->nx);

    check_nan_2d(p->temp->current, p->energy->nx ,p->energy->ny, "temperature");
    check_nan_2d(p->energy->current, p->energy->nx ,p->energy->ny, "energy");
    check_monotoic_up(p->energy->current, p->energy->nx ,p->energy->ny);
    check_monotoic_up(p->temp->current, p->energy->nx ,p->energy->ny);
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
double ***build_matrix_A(Coordinate *coor, Data *vol,Data *diff, double dt, Constants* constants, Quantity * eng) {
    int i = 0,j = 0, k;
    int i1,i0,j1,j0;
    int nxp = coor->nxp;
    int nyp = coor->nyp;
    int nx = vol->nx;
    int ny = vol->ny;
    double **volume = vol->values;
    double **R = coor->R;
    double **Z = coor->Z;
    double **D = diff->values;
    double **energy = eng->prev;
    double c_light = constants->c_light;
    double TH = constants->TH;
    double xl1, xl0;
    double sigma_boltz = constants->sigma_boltzman;
    double **RK_KL = malloc_2d(nxp, nyp);
    double **RL_KL = malloc_2d(nxp, nyp);
    double **ZK_KL = malloc_2d(nxp, nyp);
    double **ZL_KL = malloc_2d(nxp, nyp);
    
    double **rho1 = malloc_2d(nx, ny);
    double **rho2 = malloc_2d(nx, ny);
    double **rho3 = malloc_2d(nx, ny);
    double **rho4 = malloc_2d(nx, ny);

    double **S_sigma = malloc_2d(nx, ny);
    double **S_lambda = malloc_2d(nx, ny);

    double **B_sigma = malloc_2d(nx, ny);
    double **B_lambda = malloc_2d(nx, ny);

    double **jacobi = malloc_2d(nxp, nyp);
    double ***A = malloc_3d(nxp, nyp, 10);
    //#pragma omp parallel for default(shared)
    for(i = 0; i < nxp - 1; i++){
        for (j = 0; j < nyp - 1; j++) {
            RK_KL[i][j] = (R[i + 1][j + 1] - R[i][j] + R[i + 1][j] - R[i][j + 1]) / 2.0;
            RL_KL[i][j] = (R[i + 1][j + 1] - R[i][j] - R[i + 1][j] + R[i][j + 1]) / 2.0;

            ZK_KL[i][j] = (Z[i + 1][j + 1] - Z[i][j] + Z[i + 1][j] - Z[i][j + 1]) / 2.0;
            ZL_KL[i][j] = (Z[i + 1][j + 1] - Z[i][j] - Z[i + 1][j] + Z[i][j + 1]) / 2.0;
            // printf("hey %10e\t%10e\t%10e\t%10e\n",RK_KL[i][j],RL_KL[i][j],ZK_KL[i][j],ZL_KL[i][j]);
            //calclulate jacobian..
            jacobi[i][j] = RK_KL[i][j] * ZL_KL[i][j] - RL_KL[i][j] * ZK_KL[i][j];

        }
    }

    // these sizes are cell quantities, therefore we go KC, LC
    // we do not go throught the last imaginary cells, ofc.
    // #pragma omp parallel for default(shared)
    print_2d(D, nx, ny);
    for(i = 0; i < nx - 1; i++){
        for (j = 0; j < ny - 1; j++) {
            //the reason we take jacob(i,j) is because we already made sigma(i,j)
            // be the same index as jacob(i,j), only R is shifted. see (15)
            // B_sigma[i][j] = (R[i + 1][j + 1] + R[i + 1][j]) 
            //                 / ( D[i][j] * jacobi[i][j] + D[i + 1][j] * jacobi[i + 1][j]);

            // B_lambda[i][j] = (R[i + 1][j + 1] + R[i][j + 1]) 
            // / ( D[i][j] * jacobi[i][j] + D[i][j + 1] * jacobi[i][j + 1]);
            B_sigma[i][j] = (R[i + 1][j + 1] + R[i + 1][j]) 
                            / (jacobi[i][j]/D[i][j] + jacobi[i + 1][j] / D[i + 1][j]);

            B_lambda[i][j] = (R[i + 1][j + 1] + R[i][j + 1]) 
                            / (jacobi[i][j]/D[i][j] + jacobi[i][j + 1]/D[i][j + 1]);
        }
    }
    // exit(1);
    // BOUNDARY CONDITIONS
    // TOP WALL
    for (i = 0; i < nx; i++) {
        // data[i][m-1] = 0
        // B_sigma[i][ny-1] = 0;
        B_lambda[i][ny - 2] = 0;
    }

    // Bottom WALL - reflective
    // i and j = 0
    for (i = 0; i < nx; i++) {
        // data[i][0]     = data[i][1]; //no need
        B_sigma[i][0] = 0;
        B_lambda[i][0] = 0;
    }

    //RIGHT WALL - leapes
    for (j = 0; j < ny; j++) {
        // data[n - 1][j] = 0; nothing in our case
        B_sigma[nx-2][j] = 0;
    }
    
    //LEFT WALL - 190 ev DO DO DO
    i = 0;
    for (j = 0; j < ny; j++) {
        // something with TH here and finc
        i1 = i + 1;
        i0 = i;
        xl0 = 0.5 * sqrt( pow(R[i0 + 1][j] - R[i0][j],2) + pow(Z[i0 + 1][j] - Z[i0][j],2) );
        xl1 = 0.5 * sqrt( pow(R[i1 + 1][j] - R[i1][j],2) + pow(Z[i1 + 1][j] - Z[i1][j],2) );
        B_sigma[i][j] = (R[i + 1][j + 1] + R[i + 1][j]) 
                            / (jacobi[i][j] + jacobi[i + 1][j]) 
                            * 2.0 * (sigma_boltz * pow(TH, 4) - c_light * energy[i][j] / 4.0) / (energy[i+1][j] - energy[i][j] + 1E-30) * (xl1 + xl0);
        // xl1 = 0.5*sqrt( pow(R[i][j1+1] - R[i][j1], 2) + pow(Z[i][j1+1]-Z[i][j1] , 2));
        // xl0 = 0.5*sqrt( pow(R[i][j0+1] - R[i][j0], 2) + pow(Z[i][j0+1]-Z[i][j0] , 2));
        // B_lambda[i][j] = (R[i + 1][j + 1] + R[i][j + 1]) 
        //             / (jacobi[i][j] + jacobi[i][j + 1]) 
        //             * 2.0 * (sigma_boltz * pow(TH, 4) - c_light * energy[i][j] / 4.0) / (energy[i][j+1] - energy[i][j] + 1E-30) * (xl1 + xl0);
                    
    }
            
    // print_2d(B_lambda, nx, ny);
    // printf("\n\n");
    // print_2d(B_sigma, nx, ny);
    //maybe merge with the next loop for speed..
    //calculate the rhos
    double C_KL;
    double sqrt_sigma_ij,sqrt_sigma_i1j, sqrt_lambda_ij, sqrt_lambda_ij1;
    k = 0;
    int im, jm;
    for(i = 0; i < nx; i++){
        for (j = 0; j < ny; j++) {
            im = maximum(i - 1, 0);
            jm = maximum(j - 1, 0);
            sqrt_sigma_ij   = sqrt( maximum(B_sigma[i][j], 0)); 
            sqrt_sigma_i1j  = sqrt( maximum(B_sigma[im][j], 0 ) );
            sqrt_lambda_ij  = sqrt( maximum(B_lambda[i][j] , 0) );
            sqrt_lambda_ij1 = sqrt( maximum(B_lambda[i][jm] , 0) );
            C_KL            = RK_KL[i][j] * RL_KL[i][j] + ZK_KL[i][j] * ZL_KL[i][j];
            // printf("%10e\n",B_sigma[i][j]);
            rho1[i][j]      = sqrt_sigma_ij  * sqrt_lambda_ij  * C_KL; 
            rho2[i][j]      = sqrt_sigma_i1j * sqrt_lambda_ij1 * C_KL; 
            rho3[i][j]      = sqrt_sigma_i1j * sqrt_lambda_ij  * C_KL;
            rho4[i][j]      = sqrt_sigma_ij  * sqrt_lambda_ij1 * C_KL;
        }
    }
    //calculate the small sigma and small lambda
    double r1,r2,r3,r4;
    double z1,z2,z3,z4;
    for(i = 0; i < nx; i++){
        for (j = 0; j < ny ; j++) {
            r1 = RL_KL[i][j]     * RL_KL[i][j] ;
            r2 = RL_KL[i + 1][j] * RL_KL[i+1][j];
            r3 = RK_KL[i][j]     * RK_KL[i][j];
            r4 = RK_KL[i][j + 1] * RK_KL[i][j + 1];
            z1 = ZL_KL[i][j]     * ZL_KL[i][j];
            z2 = ZL_KL[i + 1][j] * ZL_KL[i + 1][j];
            z3 = ZK_KL[i][j]     * ZK_KL[i][j];
            z4 = ZK_KL[i][j + 1] * ZK_KL[i][j + 1];

            S_sigma[i][j]  = B_sigma[i][j]  * ( r1 + z1 + r2 + z2) / 2.0;
            S_lambda[i][j] = B_lambda[i][j] * ( r3 + z3 + r4 + z4) / 2.0;
        }
    }
  
    // the matrix is represnted in the following manner:
    // 6 7 8
    // 3 4 5
    // 0 1 2    
    double vol_1;
    double a5, a6, a7, a8;
    int p = 0;
    k = 0;
    for(i = 0; i < nx - 1; i++){
        //k is equivilent to i - 1. instead of MAX(i,0)...
        for (j = 0; j < ny - 1; j++) {
            im = maximum(i - 1, 0);
            jm = maximum(j - 1, 0);
            vol_1 = 1.0 / volume[i][j];
            A[i][j][4] = 1.0 - dt*vol_1*( -S_sigma[i][j] - S_sigma[im][j] - S_lambda[i][j] - S_lambda[i][jm]
                         + 0.5* (rho1[i][j] + rho2[i][j] - rho3[i][j] - rho4[i][j]) );
            a8 = 0.25 * dt * (rho3[i + 1][j] + rho4[i][j + 1]);
            a7 = dt * (S_lambda[i][j] - 0.25*(rho1[i][j] + rho2[i][j + 1] - rho3[i][j] - rho4[i][j + 1]));  
            a6 = 0.25 * dt * (rho1[im][j] + rho2[i][j + 1]);
            a5 = dt * (S_sigma[i][j] - 0.25 * (rho1[i][j] + rho2[i + 1][j] - rho3[i + 1][j] - rho4[i][j]));
            
            A[i][j][8]         = -vol_1 * a8; // equation (17d)
            A[i][j][7]         = -vol_1 * a7; // equation (17c)
            A[i][j][6]         = -vol_1 * a6; // equation (17e)
            A[i][j][5]         = -vol_1 * a5; // equation (17b)

            A[i + 1][j + 1][0] = -(1.0/volume[i + 1][j + 1]) * a8;
            A[i][j + 1][1]     = -(1.0/volume[i][j + 1])     * a7;
            A[im][j + 1][2]    = -(1.0/volume[im][j + 1])    * a6;
            A[i + 1][j][3]     = -(1.0/volume[i + 1][j])     * a5;
        }
    }

    // for(i = 0; i < nx - 1; i++){
    //     //k is equivilent to i - 1. instead of MAX(i,0)...
    //     for (j = 0; j < ny - 1; j++) {
    //         // A[i][j][4] = 1.0;
    //         // A[i][j][1] = 0.0;
    //         // A[i][j][7] = 0.0;
    //         // A[i][j][0] = 0.0;
    //         // A[i][j][8] = 0.0;
    //         // A[i][j][2] = 0.0;
    //         // A[i][j][6] = 0.0;
    //         // A[i][j][3] = 0.0;
    //         // A[i][j][5] = 0.0;
    //     }
    // }

    free_2d(RK_KL, nxp);
    free_2d(RL_KL, nxp);
    free_2d(ZK_KL, nxp);
    free_2d(ZL_KL, nxp);

    free_2d(rho1, nx);
    free_2d(rho2, nx);
    free_2d(rho3, nx);
    free_2d(rho4, nx);

    free_2d(B_sigma, nx);
    free_2d(B_lambda, nx);
    free_2d(S_sigma, nx);
    free_2d(S_lambda, nx);
    
    free_2d(jacobi, nx);
    // print_3d(A,nx,ny,10);
    check_nan_3d(A, nx,ny,10);
    return A;
}


/**
 * @brief Prepares the solution vector b, Ax=b
 * energy - is the prev energy i.e E(r, t = n), current = n + 1 
 * 
*/
double **build_b_vector(Quantity *E, Quantity *T, Data *opac,double**X, Constants *consts, double dt, int cyc) {
    int i, j;
    int nx = T->nx;
    int ny = T->ny;
    double c = consts->c_light;
    double sigma_factor = consts->sigma_factor;
    double **b = malloc_2d(nx, ny);
    double **energy = E->prev;
    double **temp = T->prev;
    double **opacity =  opac->values;
    double src = 0.0;

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            if (X[i][j] < 0.5 && dt * consts->c_light < 10) {
                src = consts->source;
            } else {
                src = 0.0;
            }
            b[i][j] = energy[i][j] + sigma_factor * opacity[i][j] * dt * c * consts->a_rad * pow(temp[i][j],4) + src * c * dt;
        }
    }
    return b;
}

void update_prev_values(Problem * p) {
    int i, j;
    int X = p->temp->nx;
    int Y = p->temp->ny;

    double **t_prev = p->temp->prev;
    double **t_curr = p->temp->current;

    double **e_prev = p->energy->prev;
    double **e_curr = p->energy->current;

    for (i = 0; i < X; i++) {
        for (j = 0; j < Y; j++) {
            t_prev[i][j] = t_curr[i][j];
            e_prev[i][j] = e_curr[i][j];
        }
    }
}