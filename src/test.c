#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>


int main (int argc, char* argv[]) {
    int i = 0,j = 0,k;
    int K_max = atoi(argv[1]);
    int L_max = atoi(argv[2]);
    
    double **rho1  = malloc_2d(K_max,L_max);
    double **rho2  = malloc_2d(K_max,L_max); 
    double **rho3  = malloc_2d(K_max,L_max);
    double **rho4  = malloc_2d(K_max,L_max);
    double **rho11  = malloc_2d(K_max,L_max);
    double **rho22  = malloc_2d(K_max,L_max); 
    double **rho33  = malloc_2d(K_max,L_max);
    double **rho44  = malloc_2d(K_max,L_max);
    double **RK_KL = malloc_2d(K_max,L_max);
    double **RL_KL = malloc_2d(K_max,L_max);
    double **sigma = malloc_2d(K_max,L_max);
    double **lambda= malloc_2d(K_max,L_max);
    //!--------- now we calculate the RHO [1-4] ---------!
    //check two schemes.
    double C_KL;
    double sqrt_sigma_ij,sqrt_sigma_i1j, sqrt_lambda_ij, sqrt_lambda_ij1;
    for(i = 0; i < K_max - 1; i++){
        for (j = 0; j < L_max - 1; j++) {
            sigma[i][j] = (i * j) % 7;
            lambda[i][j] = (i * 3*j) % 12;
            RK_KL[i][j] = ( i + j +2) %11;
            RL_KL[i][j] = ( 3*i + j +2) %5;
        }
    }
    k = 0;

    double start = omp_get_wtime();
    for(i = 0; i < K_max - 1; i++){
       // j = 0;
        //sqrt_lambda_ij1 = sqrt( lambda[i][j - 1]);
        //sqrt_lambda_ij1 = sqrt_lambda_ij  = sqrt( lambda[i][j]);
        
        for (j = 0; j < L_max - 1; j++) {
            /*if ( i != 0) {
                sqrt_sigma_i1j  = sqrt( sigma[i - 1][j] );
            } else {

            }*/
            sqrt_sigma_i1j  = sqrt( sigma[k][j] );
            sqrt_sigma_ij   = sqrt( sigma[i][j]); 
            sqrt_lambda_ij  = sqrt( lambda[i][j]);
            C_KL            = RK_KL[i][j] * RL_KL[i][j];
            rho1[i][j]      = sqrt_sigma_ij  * sqrt_lambda_ij  * C_KL; 
            rho2[i][j]      = sqrt_sigma_i1j * sqrt_lambda_ij1 * C_KL; 
            rho3[i][j]      = sqrt_sigma_i1j * sqrt_lambda_ij  * C_KL;
            rho4[i][j]      = sqrt_sigma_ij  * sqrt_lambda_ij1 * C_KL;
            sqrt_lambda_ij1 = sqrt_lambda_ij;
        }
        k = i;
    }
   
   // printf("1. Time elapsed in ms: %f\n", omp_get_wtime() - start);
    k = 0;
   // start = omp_get_wtime();
    for(i = 0; i < K_max - 1; i++){
       // j = 0;
        //sqrt_lambda_ij1 = sqrt( lambda[i][j - 1]);
        //sqrt_lambda_ij1 = sqrt_lambda_ij  = sqrt( lambda[i][j]);

        for (j = 0; j < L_max - 1; j++) {
            sqrt_sigma_ij   = sqrt( sigma[i][j]); 
            sqrt_sigma_i1j  = sqrt( sigma[k][j] );
            sqrt_lambda_ij  = sqrt( lambda[i][j]);
            sqrt_lambda_ij1 = sqrt( lambda[i][j - 1]);
            C_KL            = RK_KL[i][j] * RL_KL[i][j];
            rho11[i][j]      = sqrt_sigma_ij  * sqrt_lambda_ij  * C_KL; 
            rho22[i][j]      = sqrt_sigma_i1j * sqrt_lambda_ij1 * C_KL; 
            rho33[i][j]      = sqrt_sigma_i1j * sqrt_lambda_ij  * C_KL;
            rho44[i][j]      = sqrt_sigma_ij  * sqrt_lambda_ij1 * C_KL;
        }
        k = i;
    }
    for(i = 0; i < K_max - 1; i++){
        for (j = 0; j < L_max - 1; j++) {
            if ( rho11[i][j] != rho1[i][j] ||
            rho22[i][j] != rho2[i][j] ||
            rho33[i][j] != rho3[i][j] ||
            rho44[i][j] != rho4[i][j]) {
                printf ("NO GOOD\n");
            }
        }
    }

    return 0;
}
