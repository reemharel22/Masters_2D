#include "problem.h"
#include "utils.h"
inline double calculate_C(double R_k,double R_l) {
    return R_k * R_l;
}


/****
 * After each time step, we update the time parameters and check if we should stop
 * the calculation.
*/
int update_time(Time *t){
    if (t->time_passed == t->time_stop) {
        return 1;
    }
    t->cycle += 1;
    t->time_passed += t->dt;

    return 0;
}

void do_timestep(Problem *p) {
    int i,j;
    
    return 0.0;
}

/***
 * Building the matrix A by kershaw scheme
 * note:
 * One can alter the algorithm to save memory usage, and maybe gain more speed
 * **/
double ***build_matrix_A(Coordinate *coor, Volume *vol,Diff_Coeff *diff, double dt) {
    int i,j,k;

    int K_max = coor->K_max, L_max = coor->L_max, KC_max = vol->KC_max, LC_max = vol->LC_max;
    
    double **volume = vol->volume;
    double **R = coor->R;
    double **Z = coor->Z;
    double **D = diff->D;

    double RK_KL[K_max][L_max];
    double RL_KL[K_max][L_max];

    double rho1[K_max][L_max];
    double rho2[K_max][L_max];
    double rho3[K_max][L_max];
    double rho4[K_max][L_max];

    double sigma_B[K_max][L_max], lambda_B[K_max][L_max];
    double jacobi[K_max][L_max];
    double ***A = malloc_3d(K_max, L_max, 10);
    //now we start the fun. first gotta build determinant and RKKL RLKL
    //TODO TEST WHERE TO CALCULATE C. IN THE FIRST LOOP. OR LAST.
    double rkl,rkl1;
    //OUR RL_KL indexing starts from 0,0 as oppose to kershaw
    //we go through the bottomed imaginary cells and not topped i cells...
    //#pragma omp parallel for default(shared) private (rkl,rkl1)
    for(i = 0; i < K_max - 1; i++){
        for (j = 0; j < L_max - 1; j++) {
            rkl         = R[i + 1][j + 1] - R[i][j];
            rkl1        = R[i + 1][j] - R[i][j + 1];
            RK_KL[i][j] = (rkl + rkl1) / 2.0;
            RL_KL[i][j] = (rkl - rkl1) / 2.0;
            //calclulate jacobian..
        }
    }
    
    //#pragma omp parallel for default(shared) private (rkl,rkl1)
    for(i = 0; i < K_max; i++){
        for (j = 0; j < L_max; j++) {
            sigma_B[i][j]  = (R[i + 1][j + 1] + R[i + 1][j]) 
            / (jacobi[i + 1][j + 1]);
            lambda_B[i][j] = 5;
        }
    }
    
    return A;
}