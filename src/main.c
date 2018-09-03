#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include "problem.h"
#include "initialize.h"


int main (int argc, char* argv[]) {
    // !--------------------- TIME ----------------!
    double dt = 0, t0 = 0, time_passed = 0, time_stop = 0;
    int cycle = 0;
    // !--------------------- END TIME ----------------!

    //!------------------- VERTEX ---------------!
    // number of vertex point
    int K_max, L_max;
    //mesh coordinates
    double **R = NULL,**Z = NULL;
    //!------------------- END VERTEX -----------!

    // !------------------- CELL ----------------!
    //number of cell points, which is suppose to be vertex points -1...
    int KC_max, LC_max;
    //Energy, flux etc.. 
    double **E_current = NULL, **E_old = NULL;
    double **volume = NULL;
    // !------------------- END CELL ------------!

    double ***A = NULL;
    //!----------------- END PARAMATERS-------------------!

    init(&K_max, &L_max, &KC_max,&LC_max, &dt,&t0, &cycle,&time_passed,&time_stop,&R, &Z, &E_current, &E_old, &A,&volume);
    
    do {

    } while( !update_time(&dt,&time_passed,&time_stop,&cycle) );
    
    clean_prog(K_max, KC_max, L_max, R, Z, E_current, E_old, A,volume);
    return 0;
}
