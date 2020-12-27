#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include "problem.h"
#include "initialize.h"
#include "structs.h"
/**
 * @file main.c
 */

/*
*
* @brief Initializes through reading from a datafile.
**/

/**
 * @brief The main, executes the whole scheme.
 * Firstly we initialize all of our structs and parameters according to the input
 * after in a do while loop it executes a time step, afterwards 
 * it checks if one should stop.
 * end by cleaning structs.
 * Some things:
 * 1. We are saving the init rho although we dont really need to.
 * 2. Problem holds the whole structs and not points i.e they are all one after another, this might create a problem.
 * 3. Density can be replacet to a single small array cause it is constant
 */

int main (int argc, char* argv[]) {
    char * str;
    FILE *fenergy, *ftemp, *fx, *fy;
    Problem *p = malloc(sizeof(struct Problem));
    p->diff_coeff = malloc(sizeof(struct Data));
    p->vol = malloc(sizeof(struct Data));
    p->rho = malloc(sizeof(struct Data));
    p->opacity = malloc(sizeof(struct Data));
    p->heat_cap = malloc(sizeof(struct Data));

    p->diag = malloc(sizeof(struct Diagnostics));
    p->mats = malloc(sizeof(struct Materials));
    
    p->temp = malloc(sizeof(struct Quantity));
    p->energy = malloc(sizeof(struct Quantity));

    p->time = malloc(sizeof(struct Time));

    p->constants = malloc(sizeof(struct Constants));
    p->coor = malloc(sizeof(struct Coordinate));
    if (argc == 1) {
        str = "datafile";
    }
    else {
        str = argv;
    }
    init(p, str);
    fenergy = fopen("data/energy.txt", "w");
    ftemp = fopen("data/temperature.txt", "w");
    fx = fopen("data/x.txt", "w");
    fy = fopen("data/y.txt", "w");
            diagnostic_energy_1d(p->energy, p->constants->a_rad,p->time->time_passed, 1, fenergy);

    do {
        if (p->time->cycle % 1 == 0) {
            printf("\n########## Starting time step number: %d ########## \n",p->time->cycle);
        }

        do_timestep(p);
        // diagnostic_energy_1d(p->energy, p->constants->a_rad,p->time->time_passed, 1, fenergy);
        // diagnostic_temp_1d(p->temp,p->time->time_passed * p->constants->c_light, 1, ftemp);
        diagnostic_energy_1d(p->energy, p->constants->a_rad,p->time->time_passed*p->constants->c_light , 1, fenergy);
        diagnostic_energy_1d(p->temp, p->constants->a_rad,p->time->time_passed*p->constants->c_light , 1, ftemp);
        if (p->time->cycle % 1 == 0) {
            printf("\n@@@@@@@@@@ Done time step number: %d @@@@@@@@@@ \n",p->time->cycle);
        }
        // if (p->time->cycle == 3) 
        //     exit(1);
    } while( !update_time(p->time, p->temp) );
    
    clean_prog(p);
    return 0;
}

