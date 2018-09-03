#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include "problem.h"
#include "initialize.h"
#include "structs.h"

int main (int argc, char* argv[]) {
    Problem prob;

    init(&prob);
    
    do {
        
    } while( !update_time(&prob.time) );
    
    clean_prog(&prob);
    return 0;
}
