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
    double *L, *U, *mainD, *solve;
    nx = p->temp->nx;
    ny = p->temp->ny;
    update_prev_values(p);

    //BC FOR OLSON
    // for (j = 0; j < ny; j++){
    //     p->energy->prev[0][j] = pow(p->constants->TH, 4) * p->constants->a_rad;
    //     p->energy->current[0][j] = pow(p->constants->TH, 4) * p->constants->a_rad;
    //     p->temp->prev[0][j]   = p->constants->TH; //pow(p->constants->TH,4);// * p->constants->a_rad;
    //     p->temp->current[0][j]   = p->constants->TH; //pow(p->constants->TH,4);// * p->constants->a_rad;
    // }
    // printf("")
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

    A = build_matrix_A(p->coor, p->vol, p->diff_coeff,p->opacity, p->time, p->constants, p->energy);
    // print_3d(A, p->temp->nx, p->temp->ny, 10);
    b = build_b_vector(p->energy, p->temp, p->opacity,p->coor->R,p->vol, p->constants, p->time->dt, p->time->cycle);
    // print_2d(b,p->energy->nx ,p->energy->ny );

    if (VERBOSE)
        printf("Solving the matrix\n");

    // Now we will solve the tridiagonal
    L = malloc_1d(ny);
    mainD = malloc_1d(ny);
    U = malloc_1d(ny);
    solve = malloc_1d(ny);
    // printf("A\n");

    // for (size_t i = 0; i < ny - 1; i++)
    // {
    //     L[i]     = A[1][i][5]; 
    //     U[i]     = A[1][i][1];
    //     mainD[i] = A[1][i][3];
    //     solve[i] = b[1][i]; 

    // }
    // read_file("/home/reemh/leeor_olson/fort.80", mainD, ny);
    // printf("AAA\n");
    // read_file("/home/reemh/leeor_olson/fort.81", L, ny);
    // read_file("/home/reemh/leeor_olson/fort.82", U, ny);
    // read_file("/home/reemh/leeor_olson/fort.83", solve, ny);
    // // mat_mul(A, p->energy->current, b, 0.001, p->energy->nx, p->energy->ny);
    // for (size_t i = 0; i < ny-1; i++)
    // {
    //     // print(U[i]);
    // }
    
    // solveTriagonal(ny - 1, solve, L, U, mainD);
    // // print(solve[0]);
    // for (size_t j = 1; j < ny; j++)
    // {
        
    //     p->energy->current[1][j] = solve[j];
    // }
    
    jacobi_method_naive(1000, p->energy->nx, p->energy->ny, 0.0001, A, p->energy->current, b);


    free_3d(A,p->energy->nx + 1, p->energy->ny + 1);
    check_nan_2d(p->energy->current, p->energy->nx ,p->energy->ny, "energy");
    // print(p->energy->current[1][1]);
    if (VERBOSE)
        printf("Calculating Temperature\n");
    calculate_temperature(p->temp, p->energy, p->constants, p->opacity, p->heat_cap, p->time->dt);

    if (VERBOSE)
        printf("Updating the current prev values\n");
    check_nan_2d(p->temp->current, p->energy->nx ,p->energy->ny, "temperature");
    
    free_2d(b, p->energy->nx);

    check_nan_2d(p->temp->current, p->energy->nx ,p->energy->ny, "temperature");
    check_nan_2d(p->energy->current, p->energy->nx ,p->energy->ny, "energy");
    check_monotoic_up(p->energy->current, p->energy->nx ,p->energy->ny);
    check_monotoic_up(p->temp->current, p->energy->nx ,p->energy->ny);
    // check_bc(p->energy, p->coor,p->diff_coeff, p->constants);
    return;
    
}

/**
 * @brief Builds the matrix A by Kershaw scheme.
 * @param coor  the coordinates struct.
 * @param vol   the volume struct.
 * @param diff  the diffusion coefficient struct.
 * @param dt    the size of timestep.
 * 1. Calculates XK_KL XL_KL Z_KL Z_LK by equation 15c and 15d.
 * 2. Calculates the big SIGMA and big LAMBDA, although they are named sigma, lambda
 *    they are still big ones. (by equation 15a 15b).
 * 3. Calculates Rho[1~4] by equation 17e
 * 4. Calculates small sigma and small lambda.
 * 5. Calculates the matrix A by 17a-17f.
 */
double ***build_matrix_A(Coordinate *coor, Data *vol,Data *diff, Data*opac, Time*time , Constants* constants, Quantity * eng) {
    int i = 0,j = 0, k;
    int i1,i0,j1,j0;
    int nxp = coor->nxp;
    int nyp = coor->nyp;
    int nx = vol->nx;
    int ny = vol->ny;
    double **volume = vol->values;
    double **R = coor->R;
    double **Z = coor->Z;
    double **X = coor->X;
    double **Y = coor->Y;
    FILE* fp_sigma;
    FILE* fp_lambda;
    double **D = diff->values;
    double dt = time->dt;
    double **energy = eng->prev;
    double c_light = constants->c_light;
    double TH = constants->TH;
    double xl1, xl0;
    double **opacity = opac->values;
    double sigma_boltz = constants->sigma_boltzman;
    double **XK_KL = malloc_2d(nxp, nyp);
    double **XL_KL = malloc_2d(nxp, nyp);
    double **YK_KL = malloc_2d(nxp, nyp);
    double **YL_KL = malloc_2d(nxp, nyp);
    
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
            XK_KL[i][j] = (X[i + 1][j + 1] - X[i][j] + X[i + 1][j] - X[i][j + 1]) / 2.0;
            XL_KL[i][j] = (X[i + 1][j + 1] - X[i][j] - X[i + 1][j] + X[i][j + 1]) / 2.0;

            YK_KL[i][j] = (Y[i + 1][j + 1] - Y[i][j] + Y[i + 1][j] - Y[i][j + 1]) / 2.0;
            YL_KL[i][j] = (Y[i + 1][j + 1] - Y[i][j] - Y[i + 1][j] + Y[i][j + 1]) / 2.0;
            //calclulate jacobian..
            jacobi[i][j] = XK_KL[i][j] * YL_KL[i][j] - XL_KL[i][j] * YK_KL[i][j];
        }
    }

    
    // print_2d(jacobi, nx, ny);
    fp_lambda = fopen("data/jacobi.txt", "a");
    
    fprintf(fp_lambda, "%10e ", time->time_passed);
    for (j = 0; j < ny; j++) {
        fprintf(fp_lambda, "%10e ", jacobi[1][j]);
    }
    fprintf(fp_lambda, "\n");
    fclose(fp_lambda);
    // these sizes are cell quantities, therefore we go KC, LC
    // we do not go throught the last imaginary cells, ofc.
    // #pragma omp parallel for default(shared)
    // print_2d(D, nx, ny);
    for(i = 0; i < nx -1; i++){
        for (j = 0; j < ny-1 ; j++) {
            //the reason we take jacob(i,j) is because we already made sigma(i,j)
            // be the same index as jacob(i,j), only R is shifted. see (15)
            // B_sigma[i][j] = (R[i + 1][j + 1] + R[i + 1][j]) 
            //                 / ( D[i][j] * jacobi[i][j] + D[i + 1][j] * jacobi[i + 1][j]);

            // B_lambda[i][j] = (R[i + 1][j + 1] + R[i][j + 1]) 
            // / ( D[i][j] * jacobi[i][j] + D[i][j + 1] * jacobi[i][j + 1]);
            B_sigma[i][j] = (R[i + 1][j + 1] + R[i + 1][j]) 
                            / (jacobi[i][j] + jacobi[i + 1][j]) * D[i][j];
            // printf("sigma: %d %d %10e %10e %10e %10e\n", i,j, B_sigma[i][j], jacobi[i][j], jacobi[i + 1][j], D[i][j]);
            B_lambda[i][j] = (R[i + 1][j + 1] + R[i][j + 1]) 
                            / (jacobi[i][j] + jacobi[i][j + 1]) * D[i][j];
        }
    }

    // exit(1);
    // BOUNDARY CONDITIONS
    // TOP WALL
    for (i = 0; i < nx; i++) {
        B_lambda[i][ny - 2] = 0;
    }

    // Bottom WALL 
    // i and j = 0
    j = 0;
    for (i = 0; i < nx; i++) {
        B_sigma[i][j] = 0;
        B_lambda[i][j] = 0;
        j1 = j + 1;
        j0 = j;
        xl0 = 0.5 * sqrt( pow(X[i][j0 + 1] - X[i][j0],2) + pow(Y[i][j0 + 1] - Y[i][j0],2) );
        xl1 = 0.5 * sqrt( pow(X[i][j1 + 1] - X[i][j1],2) + pow(Y[i][j1 + 1] - Y[i][j1],2) );
        if (nx == 3) {
        B_lambda[i][j] = (R[i + 1][j + 1] + R[i][j + 1]) 
                    / (jacobi[i][j] + jacobi[i][j + 1]) 
                    * 2.0 * (sigma_boltz * pow(TH, 4) - c_light * energy[i][j+1] / 4.0) / (energy[i][j] - energy[i][j+1] ) * (xl1 + xl0);
        }
    }

    //RIGHT WALL - leapes
    for (j = 0; j < ny; j++) {
        // data[n - 1][j] = 0; nothing in our case
        // B_sigma[nx-2][j] = 0;
    }
    
    //LEFT WALL - 190 ev DO DO DO
    i = 0;
    for (j = 0; j < ny; j++) {
        // something with TH here and finc
        i1 = i + 1;
        i0 = i;
        xl0 = 0.5 * sqrt( pow(X[i0 + 1][j] - X[i0][j],2) + pow(Y[i0 + 1][j] - Y[i0][j],2) );
        xl1 = 0.5 * sqrt( pow(X[i1 + 1][j] - X[i1][j],2) + pow(Y[i1 + 1][j] - Y[i1][j],2) );
        B_sigma[i][j] = 0;
        // B_lambda[i][j] = 0;
        if (ny == 3) {
            B_sigma[i][j] = (R[i + 1][j + 1] + R[i + 1][j]) 
                            / (jacobi[i][j] + jacobi[i + 1][j]) 
                            * 2.0 * (sigma_boltz * pow(TH, 4) - c_light * energy[i+1][j] / 4.0) / (energy[i][j] - energy[i+1][j]) * (xl1 + xl0);
        }
        // B_sigma[i][j] = 122432.11547629160*122432.11547629160;
        // printf("B_SIGMA i: %d j: %d value:%10e\n",i,j,B_sigma[i][j]);
        // xl1 = 0.5*sqrt( pow(R[i][j1+1] - R[i][j1], 2) + pow(Z[i][j1+1]-Z[i][j1] , 2));
        // xl0 = 0.5*sqrt( pow(R[i][j0+1] - R[i][j0], 2) + pow(Z[i][j0+1]-Z[i][j0] , 2));
    }
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
            C_KL            = XK_KL[i][j] * XL_KL[i][j] + YK_KL[i][j] * YL_KL[i][j];
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
            r1 = XL_KL[i][j]     * XL_KL[i][j] ;
            r2 = XL_KL[i + 1][j] * XL_KL[i+1][j];
            r3 = XK_KL[i][j]     * XK_KL[i][j];
            r4 = XK_KL[i][j + 1] * XK_KL[i][j + 1];
            z1 = YL_KL[i][j]     * YL_KL[i][j];
            z2 = YL_KL[i + 1][j] * YL_KL[i + 1][j];
            z3 = YK_KL[i][j]     * YK_KL[i][j];
            z4 = YK_KL[i][j + 1] * YK_KL[i][j + 1];

            S_sigma[i][j]  = B_sigma[i][j]  * ( r1 + z1 + r2 + z2) / 2.0;
            S_lambda[i][j] = B_lambda[i][j] * ( r3 + z3 + r4 + z4) / 2.0;
        }
    }
        // S_sigma[0][1] = S_sigma[nx-2][1] = 2*S_sigma[nx-3][1] ;

    fp_sigma  = fopen("data/sigma.txt", "a");
    fp_lambda = fopen("data/lambda.txt", "a");

    fprintf(fp_sigma, "%10e ", time->time_passed);
    for (i = 0; i < nx - 1; i++) {
        fprintf(fp_sigma, "%10e ", S_sigma[i][1]);
    }
    fprintf(fp_sigma, "\n");

    fprintf(fp_lambda, "%10e ", time->time_passed);
    for (j = 0; j < ny - 1; j++) {
        fprintf(fp_lambda, "%10e ", S_lambda[1][j]);
    }
    fprintf(fp_lambda, "\n");
    fclose(fp_lambda);
    fclose(fp_sigma); 
//LEEOR:
// 1 4 7
// 2 5 8
// 3 6 9


// X X   17e
// X 17a 17c
// X 17b 17d
    // the matrix is represnted in the following manner:
    // 6 7 8
    // 3 4 5
    // 0 1 2    
    double vol_1;
    double b_17, c_17, d_17, e_17;
    int p = 0;
    double cdt = c_light * dt;
    k = 0;
    for(i = 0; i < nx - 1; i++){
        //k is equivilent to i - 1. instead of MAX(i,0)...
        for (j = 0; j < ny - 1; j++) {
            im = maximum(i - 1, 0);
            jm = maximum(j - 1, 0);
            vol_1 = 1.0 / volume[i][j];
            A[i][j][4] = (1.0 * volume[i][j] + volume[i][j]*cdt * opacity[i][j]) - dt*( -S_sigma[i][j] - S_sigma[im][j] - S_lambda[i][j] - S_lambda[i][jm]
                         + 0.5* (rho1[i][j] + rho2[i][j] - rho3[i][j] - rho4[i][j]) );

            b_17 = S_sigma[i][j]  - 0.25 * (rho1[i][j] + rho2[i + 1][j] - rho3[i + 1][j] - rho4[i][j]);
            c_17 = S_lambda[i][j] - 0.25 * (rho1[i][j] + rho2[i][j + 1] - rho3[i][j]     - rho4[i][j + 1]);  
            d_17 = -0.25 * (rho3[i + 1][j] + rho4[i][j + 1]);
            e_17 =  0.25 * (rho1[im][j] + rho2[i][j + 1]);
            
            A[i][j][1]         = -dt  * b_17;// * vol_1; // equation (17b)
            A[i][j][5]         = -dt  * c_17;// * vol_1; // equation (17c)
            A[i][j][2]         = -dt  * d_17;// * vol_1; // equation (17d)
            A[i][j][8]         = -dt  * e_17;// * vol_1; // equation (17e)

            A[i+1][j]  [7] = A[i][j][1]; //* volume[i][j]/volume[i+1][j];
            A[i]  [j+1][3] = A[i][j][5]; // * volume[i][j]/volume[i]  [j+1];
            A[i+1][j+1][6] = A[i][j][2]; //* volume[i][j]/volume[i+1][j+1];
            A[im] [j+1][0] = A[i][j][8]; //* volume[i][j]/volume[im] [j+1];

        }
    }
    for (size_t j = 0; j < ny - 1; j++)
    {
        // print(A[1][j][4]);
    }
    
    free_2d(XK_KL, nxp);
    free_2d(XL_KL, nxp);
    free_2d(YK_KL, nxp);
    free_2d(YL_KL, nxp);

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
double **build_b_vector(Quantity *E, Quantity *T, Data *opac,double**X, Data*vol, Constants *consts, double dt, int cyc) {
    int i, j;
    int nx = T->nx;
    int ny = T->ny;
    double c = consts->c_light;
    double sigma_factor = consts->sigma_factor;
    double **b = malloc_2d(nx, ny);
    double **energy = E->prev;
    double **temp = T->prev;
    double **opacity =  opac->values;
    double ** volume = vol->values;
    double src = 0.0;

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            if (X[i][j] < 0.5 && dt * consts->c_light < 10) {
                src = consts->source;
            } else {
                src = 0.0;
            }
            b[i][j] = energy[i][j] + sigma_factor * opacity[i][j] * dt * c * consts->a_rad * pow(temp[i][j], 4) + src * c * dt;
            b[i][j] = b[i][j] * volume[i][j];
        }
    }
    
    // print(b[0][1]);
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

    for (i = 1; i < X-1; i++) {
        for (j = 1; j < Y-1; j++) {
            t_prev[i][j] = t_curr[i][j];
            e_prev[i][j] = e_curr[i][j];
        }
    }
}

void check_bc(Quantity*E, Coordinate* coord, Data*diff, Constants *consts) {
    int i,j;
    double ** e = E->current;
    double ** x = coord->X;
    int nx = E->nx;
    int ny = E->ny;
    double f = 0.0;// malloc_2d(nx, ny);   
    double TH = consts->TH;
    double sigma_boltz = consts->sigma_boltzman;
    double c = consts->c_light;
    double bc_val;
    double tmp = 0.0;
    double **D = diff->values;

    bc_val =  consts->sigma_boltzman* pow(TH,4);
    i = 0;
    for (j = 0; j < ny; j++) {
        // e[i+1][j] = 138064.57662909097;
        f = -D[i][j] * (e[i+1][j] - e[i][j]) / (x[i+1][j] - x[i][j]);
        tmp = c * e[i][j] / 4.0 + f / 2.0;
        if (abs(tmp - bc_val) / abs(tmp + bc_val) > 0.01) {
            printf("BC NOT WORKING!\n");
            printf("E: %10e\nF: %10e\nMY_BC: %10e\nREAL_BC:%10e\n", e[i][j], f, tmp,bc_val);
        }
    }
}