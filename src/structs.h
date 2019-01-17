#ifndef struct_H_
#define struct_H_
// !-------------- 
//THE SIZE HERE IS NOT THE REAL SIZE BUT THE SIZE INCLUDE THE FAKE CELLS!!!!!!1
// --------------!

/***
 * @brief The diffusion coefficient struct, contains the size and the data. Cell quantity.
 * real data is 1.... to KC_max/LC_max - 1
 */
typedef struct Diff_Coeff{
    int KC_max; //!< Number of values in the K (X) direction. INCLUDES THE IMAGINARY CELLS.
    int LC_max; //!< Number of values in the L (Y) direction. INCLUDES THE IMAGINARY CELLS.
    double ** D;//!< Data values of the diffusion coefficients.
}Diff_Coeff;

/***
 * @brief The Volume struct, contains the size and the data. This is a cell quantity.
 * real data is 1.... to KC_max/LC_max - 1
 */
typedef struct Volume {
    int KC_max;//!< Number of values in the K (X) direction. INCLUDES THE IMAGINARY CELLS.
    int LC_max;//!< Number of values in the L (Y) direction. INCLUDES THE IMAGINARY CELLS.
    double **volume;//!< Data values of the volume.
}Volume;

/***
 * @brief The Energy struct, contains the size and the data. This is a cell quantity.
 * real data is 1.... to KC_max/LC_max - 1
 */
typedef struct Energy {
    int KC_max;         //!< Number of values in the K (X) direction. INCLUDES THE IMAGINARY CELLS.
    int LC_max;         //!< Number of values in the L (Y) direction. INCLUDES THE IMAGINARY CELLS.
    double **E_current; //!< Current energy values.
    double **E_old;     //!< Last time step energy values.
}Energy;

/***
 * @brief The Coordinates struct, contains the size and the data. This is a vertex quantity.
 * real data is 1.... to K_max/L_max -1
 */
typedef struct Coordinate {
    int K_max;  //!< Number of values in the K (X) direction. INCLUDES THE IMAGINARY CELLS.
    int L_max;  //!< Number of values in the L (Y) direction. INCLUDES THE IMAGINARY CELLS.
    double **R;//!< R axis values data.
    double **Z;//!< Z axis values data.
}Coordinate;

/***
 * @brief The Time struct, contains dt, t0 time passed, time to stop and number of cycle.
 * 
 */
typedef struct Time {
    double dt;//!< Size of time step.
    double t0;//!< Initial time.
    double time_passed;//!< time passed until this cycle.
    double time_stop;//!< time should stop the calculation.
    int cycle;       //!< number of cycle
}Time;

/***
 * @brief The Problem struct, contains all of the other structs.
 * 
 */
typedef struct Problem {
    struct Diff_Coeff diff_coeff;
    struct Volume vol;
    struct Energy eng;
    struct Coordinate coor;
    struct Time time;
}Problem;

/***
 * @brief The Constants defined in the problem
 * 
 */
typedef struct Constants {
    double a_rad;
    double pi;
    double sigma_boltzman;
    double c_light;
    
}Constants;

#endif