#ifndef calculations_H_
#define calculations_H_
#include "structs.h"

void calculate_temperature(Quantity *T, Quantity *E, Constants *consts, Data *opac, Data *heat_cap, double dt);
void calculate_opacity(Data *opacity,Data *rho, Quantity *T, Constants * consts);
#endif
