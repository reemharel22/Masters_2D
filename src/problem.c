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
 * THE TEMPERATURE ISNT REALLY T BUT V...
 *  
 */
void do_timestep(Problem *p) {
    double ***A, **b;
    apply_boundary(p->temp->prev, p->boundary_type, p->energy->nx, p->energy->ny);

    // @@@ DEBUG
    check_nan_2d(p->temp->prev,  p->energy->nx, p->energy->ny, "Temp prev");
    check_nan_2d(p->energy->prev,  p->energy->nx, p->energy->ny, "Energy prev");
    // @@@ DEBUG

    calculate_opacity(p->opacity, p->rho, p->temp, p->mats); //with prev
    // @@@ DEBUG
        check_nan_2d(p->opacity->values,  p->opacity->nx, p->opacity->ny,"Opacity");
        check_nan_2d(p->rho->values,  p->rho->nx, p->rho->ny, "Rho");
    // @@@ DEBUG

    calculate_diffusion_coefficient(p->diff_coeff, p->opacity, p->constants, 0);
    // @@@ DEBUG
        check_nan_2d(p->diff_coeff->values,  p->diff_coeff->nx, p->diff_coeff->ny, "Diff coeff");
    // @@@ DEBUG



    calculate_heat_capacity(p->heat_cap, p->rho, p->temp, p->mats);
        // @@@ DEBUG
        check_nan_2d(p->heat_cap->values,  p->heat_cap->nx, p->heat_cap->ny,"Heat capacity");
    // @@@ DEBUG

    //apply boundary on the energy and temperature
    
    //Calculating the energy.
    A = build_matrix_A(p->coor, p->vol, p->diff_coeff, p->time->dt);
    print_3d(A, p->temp->nx, p->temp->ny, 10);
    b = build_b_vector(p->energy, p->temp, p->opacity, p->constants, p->time->dt);

    exit(1);
    jacobi_method_naive(1000, p->energy->nx, p->energy->ny, 1e-10, A, p->energy->current, b);
    //print_2d(p->energy->current, p->temp->nx, p->temp->ny);
    
    free_3d(A,p->energy->nx, p->energy->ny);
    calculate_temperature(p->temp, p->energy, p->constants, p->opacity, p->heat_cap, p->time->dt);

    update_prev_values(p);
    free_2d(b, p->energy->ny);
    if (p->time->cycle == 2) {
        
    }
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
    int i = 0,j = 0, k;

    int nxp = coor->nxp;
    int nyp = coor->nyp;
    int nx = vol->nx;
    int ny = vol->ny;
    double **volume = vol->values;
    double **R = coor->R;
    double **Z = coor->Z;
    double **D = diff->values;

    double **RK_KL = malloc_2d(nxp, nyp);
    double **RL_KL = malloc_2d(nxp, nyp);
    double **ZK_KL = malloc_2d(nxp, nyp);
    double **ZL_KL = malloc_2d(nxp, nyp);
    
    double **rho1 = malloc_2d(nxp, nyp);
    double **rho2 = malloc_2d(nxp, nyp);
    double **rho3 = malloc_2d(nxp, nyp);
    double **rho4 = malloc_2d(nxp, nyp);

    double **S_sigma = malloc_2d(nxp, nyp);
    double **S_lambda = malloc_2d(nxp, nyp);

    double **B_sigma = malloc_2d(nx, ny);
    double **B_lambda = malloc_2d(nx, ny);

    double **jacobi = malloc_2d(nxp, nyp);
    double ***A = malloc_3d(nx, ny, 10);
    
    //#pragma omp parallel for default(shared)
    for(i = 0; i < nxp - 1; i++){
        for (j = 0; j < nyp - 1; j++) {
            RK_KL[i][j] = (R[i + 1][j + 1] - R[i][j] + R[i + 1][j] - R[i][j + 1]) / 2.0;
            RL_KL[i][j] = (R[i + 1][j + 1] - R[i][j] - R[i + 1][j] + R[i][j + 1]) / 2.0;

            ZK_KL[i][j] = (Z[i + 1][j + 1] - Z[i][j] + Z[i + 1][j] - Z[i][j + 1]) / 2.0;
            ZL_KL[i][j] = (Z[i + 1][j + 1] - Z[i][j] - Z[i + 1][j] + Z[i][j + 1]) / 2.0;
            
            //calclulate jacobian..
            jacobi[i][j] = RK_KL[i][j] * ZL_KL[i][j] - RL_KL[i][j] * ZK_KL[i][j];

        }
    }
    check_nan_2d(jacobi, nxp, nyp,"jacobi");
    check_nan_2d(RK_KL, nxp, nyp, "RKKL");
    check_nan_2d(RL_KL, nxp, nyp,"RLKL");
    check_nan_2d(ZK_KL, nxp, nyp,"ZKKL");   
    check_nan_2d(ZL_KL, nxp, nyp,"ZLKL");

    // these sizes are cell quantities, therefore we go KC, LC
    // we do not go throught the last imaginary cells, ofc.
    // #pragma omp parallel for default(shared)
    for(i = 0; i < nx - 1; i++){
        for (j = 0; j < ny - 1; j++) {
            //the reason we take jacob(i,j) is because we already made sigma(i,j)
            // be the same index as jacob(i,j), only R is shifted. see (15)
            B_sigma[i][j] = (R[i + 1][j + 1] + R[i + 1][j]) 
            / ( D[i][j] * jacobi[i][j] + D[i + 1][j] * jacobi[i + 1][j]);

            B_lambda[i][j] = (R[i + 1][j + 1] + R[i][j + 1]) 
            / ( D[i][j] * jacobi[i][j] + D[i][j + 1] * jacobi[i][j + 1]);
        }
    }
  check_nan_2d(jacobi, nxp, nyp,"B_lambda");
    //for boundary condition the sigma[0][j=0..L] = 0 and for reflective sigma[nxp][j=0..L] = 0
    for (j = 0; j < ny; j++) {
        B_sigma[nx - 1][j] = B_sigma[nx - 2][j] = B_sigma[1][j] = B_sigma[0][j] = 0;
    }
        check_nan_2d(B_sigma, nxp, nyp, "B_sigma");

    //maybe merge with the next loop for speed..
    //calculate the rhos
    double C_KL;
    double sqrt_sigma_ij,sqrt_sigma_i1j, sqrt_lambda_ij, sqrt_lambda_ij1;
    k = 0;
    int im, jm;
    for(i = 0; i < nxp - 1; i++){
        for (j = 0; j < nyp - 1; j++) {
            im = maximum(i - 1, 0);
            jm = maximum(j - 1, 0);
            sqrt_sigma_ij   = sqrt( B_sigma[i][j] ); 
            sqrt_sigma_i1j  = sqrt( B_sigma[im][j] );
            sqrt_lambda_ij  = sqrt( B_lambda[i][j] );
            sqrt_lambda_ij1 = sqrt( B_lambda[i][jm] );
            C_KL            = RK_KL[i][j] * RL_KL[i][j] + ZK_KL[i][j] * ZL_KL[i][j];
            rho1[i][j]      = sqrt_sigma_ij  * sqrt_lambda_ij  * C_KL; 
            rho2[i][j]      = sqrt_sigma_i1j * sqrt_lambda_ij1 * C_KL; 
            rho3[i][j]      = sqrt_sigma_i1j * sqrt_lambda_ij  * C_KL;
            rho4[i][j]      = sqrt_sigma_ij  * sqrt_lambda_ij1 * C_KL;
        }
    }
            check_nan_2d(rho1, nxp, nyp, "rho1");
            check_nan_2d(rho2, nxp, nyp, "rho2");
            check_nan_2d(rho3, nxp, nyp, "rho3");
            check_nan_2d(rho4, nxp, nyp, "rho4");

    //calculate the small sigma and small lambda
    double r1,r2,r3,r4;
    double z1,z2,z3,z4;
    for(i = 0; i < nxp - 1; i++){
        for (j = 0; j < nyp - 1; j++) {
            r1 = RL_KL[i][j]     * RL_KL[i][j] ;
            r2 = RL_KL[i + 1][j] * RL_KL[i+1][j];
            r3 = RK_KL[i][j]     * RK_KL[i][j];
            r4 = RK_KL[i][j + 1] * RK_KL[i][j + 1];
            z1 = ZL_KL[i][j]     * ZL_KL[i][j];
            z2 = ZL_KL[i + 1][j] * ZL_KL[i + 1][j];
            z3 = ZK_KL[i][j]     * ZK_KL[i][j];
            z4 = ZK_KL[i][j + 1] * ZK_KL[i][j + 1];

            S_sigma[i][j] = B_sigma[i][j] * ( r1 + z1  + r2 + z2)  / 2.0;
            S_lambda[i][j] = B_lambda[i][j] * ( r3 + z3 + r4 + z4) / 2.0;
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
            a6 = 0.25 *dt * (rho1[im][j] + rho2[i][j + 1]);
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

    free_2d(RK_KL, nxp);
    free_2d(RL_KL, nxp);
    free_2d(ZK_KL, nxp);
    free_2d(ZL_KL, nxp);

    free_2d(rho1, nxp);
    free_2d(rho2, nxp);
    free_2d(rho3, nxp);
    free_2d(rho4, nxp);

    free_2d(B_sigma, nxp);
    free_2d(B_lambda, nxp);
    free_2d(S_sigma, nx);
    free_2d(S_lambda, nx);
    
    free_2d(jacobi, nx);
    printf("Done with the buildA matrix\n");
    print_3d(nxp, nyp, 10, A);
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
    //    data[i][m - 1]     = data[i][m - 2]; no need
    }

    // Bottom WALL - reflective
    // i and j = 0
    for (i = 0; i < n; i++) {
        //data[i][0]     = data[i][1]; no need
    }

    //RIGHT WALL - leapes
    for (j = 0; j < m; j++) {
        //data[n - 1][j] = 0;
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
    int X = T->nx;
    int Y = T->ny;
    double c = consts->c_light;
    double sigma_factor = consts->sigma_factor;
    double **b = malloc_2d(X, Y);
    double **energy = E->prev;
    double **temp = T->prev;
    double **opacity =  opac->values;
        
    for (i = 1; i < X - 1; i++) {
        for (j = 1; j < Y - 1; j++) {
            b[i][j] = energy[i][j] + sigma_factor * opacity[i][j] * dt * c * temp[i][j];
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