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
    
    init(p, argv);
    do {
        do_timestep(p);
    } while( !update_time(p->time, p->temp) );
    
    clean_prog(p);
    return 0;
}

