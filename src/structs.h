#ifndef struct_H_
#define struct_H_
// !-------------- 
//THE SIZE HERE IS NOT THE REAL SIZE BUT THE SIZE INCLUDE THE FAKE CELLS!!!!!!1
// --------------!


/***
 * @brief the data struct which holds the 2d array that holds the values..
 */
typedef struct Data {
    int nx;
    int ny;
    double **values;
} Data;


/***
 * @brief The Energy struct, contains the size and the data. This is a cell quantity.
 * real data is 1.... to KC_max/LC_max - 1
 */
typedef struct Quantity {
    int nx;         //!< Number of values in the K (X) direction. INCLUDES THE IMAGINARY CELLS.
    int ny;         //!< Number of values in the L (Y) direction. INCLUDES THE IMAGINARY CELLS.
    double **current; //!< Current energy values.
    double **prev;     //!< Last time step energy values.
} Quantity;

/***
 * @brief The Coordinates struct, contains the size and the data. This is a vertex quantity.
 * real data is 1.... to K_max/L_max -1
 */
typedef struct Coordinate {
    int nxp;  //!< Number of values in the K (X) direction. INCLUDES THE IMAGINARY CELLS.
    int nyp;  //!< Number of values in the L (Y) direction. INCLUDES THE IMAGINARY CELLS.
    double **R;//!< R axis values data.
    double **Z;//!< Z axis values data.
} Coordinate;

/***
 * @brief The Time struct, contains dt, t0 time passed, time to stop and number of cycle.
 * 
 */
typedef struct Time {
    double dt;//!< Size of time step.
    double dt_max; // !< Max Dt.
    double t0;//!< Initial time.
    double time_passed;//!< time passed until this cycle.
    double time_stop;//!< time should stop the calculation.
    int cycle;       //!< number of cycle
    double dt_factor; // !< Dt_factor, related to update time
} Time;


/***
 * @brief The Constants defined in the problem in cgs !
 * 
 */
typedef struct Constants {
    double a_rad; //!< a radiation constant
    double c_light; //!< Speed of light
    double sigma_boltzman; //!< Sigma boltzmann
    double sigma_factor;
    double source;
    double dr;
    double dz;
    double T0;//!< Initial temperature.
} Constants;

/**
 * @brief The Material structs, holds the relevant parameters to the material.
 * The only place where this is relevant is the Cv and Opacity calculation.
*/
typedef struct Material {
    double alpha; // !< ALpha, related to kappa rossland
    double lambda; // !< Lambda related to cv
    double beta; // !< BEta related to rossaland
    double mu; // !< Mu related to rossland
    double g;  //!< g Related to rossland
    double f; // !< f related to rossland
    double init_rho; //!<density of the material
    int i_start; // !< the i index where this material starts
    int i_end;   // !< the i index where this material ends.
    int j_start;//!< the j index where this material starts.
    int j_end; //!< the j index where this material ends.
} Material;

/**
 * @brief The Material structs, holds the relevant parameters to the material.
 * The only place where this is relevant is the Cv and Opacity calculation.
*/
typedef struct Materials {
    int num_mats; //!< How many materials.
    int mat_type;
    Material *mat; //!< Material array.
} Materials;
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
    struct Data* diff_coeff;
    struct Materials* mats; 
    struct Data *vol;
    struct Data *rho;
    struct Quantity *energy;
    struct Quantity *temp;
    struct Coordinate *coor;
    struct Constants *constants;
    struct Diagnostics* diag;
    struct Data *opacity;
    struct Data *heat_cap;
    struct Time *time;
    int boundary_type;
}Problem;

#endif
