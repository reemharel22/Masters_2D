#include "utils.h"
#include "initialize.h"

/**
 * initializes the problem by parsing the file.
 * sets the K_max, L_max, aswell as their corresponding Cell.. 
 * creats the mesh, energy etc.
 * HOWEVER THIS IS NOT WHERE WE PUT THE BOUNDARY CONDITION
 * ***/
void init(int *K_max, int *L_max, int *KC_max,int *LC_max, 
double *dt, double *t0, int *cycle, double *time_passed,double *time_stop,
double ***R,double ***Z, double ***E_current, double ***E_old, double****A
 , double ***volume) {
    int i,j;
    *K_max = 78;
    *L_max = 76;
    *KC_max = *K_max - 1;
    *LC_max = *L_max - 1;
    *dt = 0.01;
    *t0 = 0;
    *time_passed = *t0;
    *time_stop = 0.02;
    *cycle = 0;
    // mallocs
    *R         = malloc_2d( *K_max , *L_max );
    *Z         = malloc_2d( *K_max , *L_max );
    *E_current = malloc_2d( *KC_max, *LC_max);
    *E_old     = malloc_2d( *KC_max, *LC_max);
    *volume    = malloc_2d( *KC_max , *LC_max );
    *A         = malloc_3d( *K_max , *L_max, 10);
    // end mallocs


    // initialize values
    init_mesh_Kershaw1(*K_max,*L_max,*R,*Z);
    for ( i = 0 ; i < *KC_max; i++) {
        for (j = 0 ; j < *LC_max; j++) {
            (*E_current)[i][j] = 0;
        }
    }
    mesh_square_volume(*volume, *R,*Z,*K_max,*L_max);
}


void init_mesh_Kershaw1(int K_max, int L_max, double **R, double **Z) {
    int i = 0, j = 0;
    for ( i = 0 ; i < K_max; i++) {
        for (j = 0 ; j < L_max; j++) {
            R[i][j] = (double) j / (L_max - 1);
        }
    }
    for ( i = 0 ; i < K_max; i++) {
        for (j = 0 ; j < L_max; j++) {
            Z[i][j] = (double) i / (K_max - 1);
        }
    }
}

void clean_prog(int K_max, int KC_max, int L_max, double **R, double **Z,
double **E_current, double **E_old, double***A,double**volume) {
    free_2d(R,K_max);
    free_2d(Z,K_max);
    free_2d(E_current,KC_max);
    free_2d(E_old,KC_max);
    free_2d(volume,KC_max);
    free_3d(A,K_max,L_max);
}
