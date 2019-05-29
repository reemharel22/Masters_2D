#include "diagnostics.h"
#include "stdio.h"


void diagnostics(Problem *p) {
    Diagnostics *d = &(p->diag);
    if (p->time.time_passed > d->time_print) {
        d->time_print *= 2;
        
    } else {

    }
    
}

void diagnostic_Quantity_2d(Quantity *Q,double t,  char*f_name) {
    int X = Q->KC_max;
    int Y = Q->LC_max;
    int i, j;
   // FILE*fp = fopen(f_name, "w");
    for (i = 0; i < X; i++) {
        for (j = 0; j < Y; j++) {
            // TODO
        }
    }
}

void diagnostic_Quantity_1d(Quantity *Q, double t, int y, char*f_name) {
    int X = Q->KC_max;
    int i;
    FILE*fp = fopen(f_name, "w");
    fprintf(fp, "%10e ", t);
    for (i = 0; i < X; i++) {
        fprintf(fp, "%10e ", Q->current[i][y]);
    }
}

void diagnostics_initial(Problem *p) {
    int i = 0;
    char* f_name = "data/Initial.txt";
    FILE*fp = fopen(f_name, "w");
    fprintf(fp, "Time Segment: {\n");
    fprintf(fp, "\tInitial time: %10e.\n\tStarting dt: %10e.\n\tdt max: %10e\n", p->time.t0, p->time.dt, p->time.dt_max);
    fprintf(fp, "\tdt factor: %10e.\n\tTime finish: %10e.\n}\n", p->time.dt_factor, p->time.time_stop);

    fprintf(fp, "Materials Segment: {\n");
    fprintf(fp, "\tNumber of Materials: %d\n", p->mats.num_mats);
    for (i = 0; i < p->mats.num_mats; i++) {
        fprintf(fp, "\tMaterial number: %d\n",i);
        fprintf(fp, "\t\tAlpha: %10e.\n",p->mats.mat[i].alpha);
        fprintf(fp, "\t\tBeta: %10e.\n",p->mats.mat[i].beta);
        fprintf(fp, "\t\tmu: %10e.\n",p->mats.mat[i].mu);
        fprintf(fp, "\t\tlambda: %10e.\n",p->mats.mat[i].lambda);
        fprintf(fp, "\t\tf: %10e.\n",p->mats.mat[i].f);
        fprintf(fp, "\t\tg: %10e.\n}",p->mats.mat[i].g);
    }

    fclose(fp);
}