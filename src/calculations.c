#include "calculations.h"
#include "math.h"

/**
 * Calculates the temperature of the material
 * 
*/
void calculate_temperature(Quantity *T, Quantity *E, Constants *consts, Data *opac, Data *heat_cap, double dt) {
    int i, j;
    int X = T->KC_max;
    int Y = T->LC_max;
    double arad, c, coeff;
    double **temp_curr, **temp_prev, **energy, **opacity, **cv;
    temp_curr = T->current;
    temp_prev = T->prev;
    energy = E->current;
    opacity = opac->values;
    cv = heat_cap->values;
    arad = consts->a_rad;
    c = consts->c_light;

    for (i = 1; i < X; i++) {
        for (j = 1; j < Y; j++) {
            coeff = (opacity[i][j] * 4.0 * pow(temp_prev[i][j], 3) * arad) / (cv[i][j]);
            temp_curr[i][j] = (temp_prev[i][j] + coeff * c * dt * energy[i][j]) 
                                / (1.0 + dt * c * coeff);
        }
    }
    return;
}

/**
 * Calculates the Opacity by rossland.
 * sigma = kappa * rho
 * 
*/
void calculate_opacity(Data *opacity,Data *rho, Quantity *T, Constants * consts) {
    int i, j;
    int X = T->KC_max;
    int Y = T->LC_max;
    double g = consts->g;
    double alpha = consts->alpha;
    double lambda = consts->lambda;
    double **temp = T->current;
    double **density = rho->values;
    double **opac = opacity->values;

    for (i = 0; i < X; i++) {
        for (j = 0; j < Y; j++) {
            opac[i][j] =  pow(density[i][j], lambda + 2.0) / (g * pow(temp[i][j], alpha));
        }
    }
}