#ifndef struct_H_
#define struct_H_
// !-------------- 
//THE SIZE HERE IS NOT THE REAL SIZE BUT THE SIZE INCLUDE THE FAKE CELLS!!!!!!1
// --------------!


/***
 * @brief the data struct which holds the 2d array that holds the values..
 */
typedef struct Data {
    int KC_max;
    int LC_max;
    double **values;
} Data;


/***
 * @brief The Energy struct, contains the size and the data. This is a cell quantity.
 * real data is 1.... to KC_max/LC_max - 1
 */
typedef struct Quantity {
    int KC_max;         //!< Number of values in the K (X) direction. INCLUDES THE IMAGINARY CELLS.
    int LC_max;         //!< Number of values in the L (Y) direction. INCLUDES THE IMAGINARY CELLS.
    double **current; //!< Current energy values.
    double **prev;     //!< Last time step energy values.
} Quantity;

/***
 * @brief The Coordinates struct, contains the size and the data. This is a vertex quantity.
 * real data is 1.... to K_max/L_max -1
 */
typedef struct Coordinate {
    int K_max;  //!< Number of values in the K (X) direction. INCLUDES THE IMAGINARY CELLS.
    int L_max;  //!< Number of values in the L (Y) direction. INCLUDES THE IMAGINARY CELLS.
    double **R;//!< R axis values data.
    double **Z;//!< Z axis values data.
} Coordinate;

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
 * @brief The Constants defined in the problem
 * 
 */
typedef struct Constants {
    double a_rad;
    double pi;
    double sigma_boltzman;
    double c_light;
    
} Constants;

/***
 * @brief The diagnostics parameters defined in the problem
 * 
 */
typedef struct Diagnostics {
    double time_print;
    char *problem_name;
    
}Diagnostics;


/***
 * @brief The Problem struct, contains all of the other structs.
 * 
 */
typedef struct Problem {
    struct Diff_Coeff diff_coeff;
    struct Volume vol;
    struct Quantity energy;
    struct Quantity temp;
    struct Coordinate coor;
    struct Diagnostics diag;
    struct Data opacity;
    struct Data heat_cap;
    struct Time time;
}Problem;

#endif