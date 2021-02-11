#include "diagnostics.h"
#include "stdio.h"


void diagnostics(Problem *p) {
    Diagnostics *d = &(p->diag);
    
}

void diagnostic_Quantity_2d(Quantity *Q,double t,  char*f_name) {
    int X = Q->nx;
    int Y = Q->ny;
    int i, j;
   // FILE*fp = fopen(f_name, "w");
    for (i = 1; i < X - 1; i++) {
        for (j = 1; j < Y - 1; j++) {

        }
    }
}

void diagnostic_position(double ** x, int nx, int ny,  FILE*fp) {
    int i, j;
    for (i = 1; i < nx - 1; i++) {
            fprintf(fp, "%10e ", x[i][1]);
    }
        fprintf(fp, "\n");

}
void diagnostic_Quantity_1d(Quantity *Q, double time, int y, FILE*fp) {
    int X = Q->nx;
    int i;
    fprintf(fp, "%10e ", time);
    for (i = 1; i < X - 1; i++) {
        fprintf(fp, "%10e ", Q->current[i][y]);
    }
    fprintf(fp, "\n");
}

void diagnostic_energy_1d(Quantity *Q, double arad, double time, int y, FILE*fp) {
    int X = Q->nx;
    int i;
    fprintf(fp, "%10e ", time);
    for (i = 1; i < X - 1; i++) {
        fprintf(fp, "%10e ", pow(Q->current[i][y]/arad, 0.25)/1160500);
    }
    fprintf(fp, "\n");
}

void diagnostic_temp_1d(Quantity *Q, double arad, double time, int y, FILE*fp) {
    int X = Q->nx;
    int i;
    fprintf(fp, "%10e ", time);
    for (i = 1; i < X - 1; i++) {
        fprintf(fp, "%10e ", arad * pow(Q->current[i][y], 4));
    }
    fprintf(fp, "\n");
}

void diagnostics_initial(Problem *p) {
    int i = 0;
    char* f_name = "data/Initial.txt";
    FILE*fp = fopen(f_name, "w");
    fprintf(fp, "Time Segment: {\n");
    fprintf(fp, "\tInitial time: %10e.\n\tStarting dt: %10e.\n\tdt max: %10e\n", p->time->t0,
     p->time->dt, p->time->dt_max);
    fprintf(fp, "\tdt factor: %10e.\n\tTime finish: %10e.\n}\n", p->time->dt_factor, p->time->time_stop);

    fprintf(fp, "Materials Segment: {\n");
    fprintf(fp, "\tNumber of Materials: %d\n", p->mats->num_mats);
    for (i = 0; i < p->mats->num_mats; i++) {
        fprintf(fp, "\tMaterial number: %d\n",i);
        fprintf(fp, "\t\tAlpha: %10e.\n",p->mats->mat[i].alpha);
        fprintf(fp, "\t\tBeta: %10e.\n",p->mats->mat[i].beta);
        fprintf(fp, "\t\tmu: %10e.\n",p->mats->mat[i].mu);
        fprintf(fp, "\t\tlambda: %10e.\n",p->mats->mat[i].lambda);
        fprintf(fp, "\t\tf: %10e.\n",p->mats->mat[i].f);
        fprintf(fp, "\t\tg: %10e.\n}",p->mats->mat[i].g);
    }

    fclose(fp);
}