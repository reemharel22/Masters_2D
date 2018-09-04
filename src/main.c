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
 */

int main (int argc, char* argv[]) {
    Problem prob;

    init(&prob);
    
    do {
        
    } while( !update_time(&prob.time) );
    
    clean_prog(&prob);
    return 0;
}
