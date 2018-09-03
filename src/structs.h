#ifndef struct_H_
#define struct_H_
// !-------------- 
//THE SIZE HERE IS NOT THE REAL SIZE BUT THE SIZE INCLUDE THE FAKE CELLS!!!!!!1
// --------------!

typedef struct Diff_Coeff{
    int KC_max, LC_max;
    double ** D;
}Diff_Coeff;

typedef struct Volume {
    int KC_max, LC_max;
    double **volume;
}Volume;

typedef struct Energy {
    int KC_max, LC_max;
    //Energy, flux etc.. 
    double **E_current;
    double **E_old;
}Energy;

typedef struct Coordinate {
    int K_max;
    int L_max;
    double **R;
    double **Z;
}Coordinate;

typedef struct Time {
    double dt;
    double t0;
    double time_passed;
    double time_stop;
    int cycle;
}Time;

typedef struct Problem {
    struct Volume vol;
    struct Energy eng;
    struct Coordinate coor;
    struct Time time;
}Problem;

#endif