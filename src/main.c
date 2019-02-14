#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include "problem.h"
#include "initialize.h"
#include "structs.h"
/**
 * @file main.c
 */

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
    Problem prob;

    init(&prob);
    do {
        
    } while( !update_time(&prob.time, &prob.temp) );
    
    clean_prog(&prob);
    return 0;
}
