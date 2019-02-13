#include "diagnostics.h"
#include "stdio.h"


void diagnostics(Problem *p) {
    Diagnostics *d = &(p->diag);
    if (p->time.time_passed > d->time_print) {
        d->time_print *= 2;
        
    } else {

    }
    
}

void diagnostic_Quantity(double t, Quantity *Q, char*f_name) {
    int X = Q->KC_max;
    int Y = Q->LC_max;
    int i, j;
    FILE*fp = fopen(f_name, "w");
    for (i = 0; i < X; i++) {
        for (j = 0; j < Y; j++) {
            // TODO
        }
    }
}