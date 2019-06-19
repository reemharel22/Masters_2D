#include "calculations.h"
#include "math.h"
#include "structs.h"
#include "utils.h"

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
            temp_curr[i][j] = pow(temp_curr[i][j] / consts->a_rad, 0.25);
        }
    }
    return;
}

/**
 * @brief Calculates the Opacity by rossland.
 * sigma = kappa * rho
 * 
*/
void calculate_opacity(Data *opacity,Data *rho, Quantity *T, Materials * mats) {
    int i, j, k = 0;
    double **temp = T->prev;
    double **density = rho->values;
    double **opac = opacity->values;
    int i_start, i_end, j_start, j_end;
    double g,lambda,alpha;
    
    for (k = 0; k < mats->num_mats; k++) {
        g = mats->mat[k].g;
        alpha = mats->mat[k].alpha;
        lambda = mats->mat[k].lambda;
        i_start = mats->mat[k].i_start;
        i_end = mats->mat[k].i_end;
        j_start = mats->mat[k].j_start;
        j_end = mats->mat[k].j_end;
        for (i = i_start; i < i_end; i++) {
            for (j = j_start; j < j_end; j++) {
                //why did i write lambda + 2
                opac[i][j] =  pow(density[i][j], lambda + 1.0) / (g * pow(temp[i][j], alpha));
            }
        }
    }
}

/**
 * @brief Calculates the Heat capacity by taking the derivative of the internal energy in respect of temperature in const volume.
 * 
*/
void calculate_heat_capacity(Data *heat_cap,Data *rho, Quantity *T, Materials * mats) {
    int i, j, k = 0;
    int X = T->KC_max;
    int Y = T->LC_max;
    double **temp = T->prev;
    double **density = rho->values;
    double **cv = heat_cap->values;
    int i_start, i_end, j_start, j_end;
    double f,beta, mu;
    double rho_mu;
    for (k = 0; k < mats->num_mats; k++) {
        f = mats->mat[k].f;
        beta = mats->mat[k].beta;
        mu = mats->mat[k].mu;
        i_start = mats->mat[k].i_start;
        i_end = mats->mat[k].i_end;
        j_start = mats->mat[k].j_start;
        j_end = mats->mat[k].j_end;
        rho_mu = pow (mats->mat[k].init_rho, 1.0 - mu);
        for (i = i_start; i < i_end; i++) {
            for (j = j_start; j < j_end; j++) {
                cv[i][j] = beta * f * pow(temp[i][j], beta - 1) * rho_mu;
            }
        }
    }
}

/**
 * 
 * Calculates the diffusion coefficent.
 * 
 * 
 * **/
void calculate_diffusion_coefficient(Data* diff, Data *opac, Constants *consts, int type) {
    int i, j;
    double c = consts->c_light;
    double **opacity = opac->values;
    double **diff_c = diff->values;
    int X = diff->KC_max;
    int Y = diff->LC_max;
    //including boundary cells
    for (i = 0; i < X; i++) {
        for (j = 0; j < Y; j++) {
            diff_c[i][j] = c / ( 3.0 * opacity[i][j]);
        }
    }
}
