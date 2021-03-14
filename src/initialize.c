#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "initialize.h"
//#include "python2.7/Python.h"
#include "diagnostics.h"
#define OLSON 1
#define SUOLSON 2
//#include "datafile.c"
/**
 * @file initialize.c
 */


/*
*
* @brief Initializes through reading from a datafile.
**/
void init_datafile(Problem *p, char* f_name) {
    FILE *fp;
    char * line = NULL;
    fp = fopen(f_name, "r");
    int num_type, num_mats;
    size_t len = 0;
    ssize_t read;
    double xmax = 0;
    double ymax = 0;
    p->constants->source = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        // set up the constants of the problem
        if(strstr(line, "sigma_factor") != NULL)
            p->constants->sigma_factor = double_reader(line, len);
       
        else if (strstr(line, "nxp") != NULL) {
            p->coor->nxp = int_reader(line, len);
        } else if(strstr(line, "nyp") != NULL) {
            p->coor->nyp = int_reader(line, len);
        } else if(strstr(line, "nx") != NULL) {
            p->vol->nx = int_reader(line, len);
        } else if(strstr(line, "ny") != NULL) {
            p->vol->ny = int_reader(line, len);
        } else if(strstr(line, "dt_0") != NULL) {
            p->time->dt = double_reader(line, len);
        } else if(strstr(line, "dt_max") != NULL) {
            p->time->dt_max = double_reader(line, len);
        } else if(strstr(line, "dt_factor") != NULL) {
            p->time->dt_factor = double_reader(line, len);
        } else if(strstr(line, "t0") != NULL) {
            p->time->t0 = double_reader(line, len);
        } else if(strstr(line, "time_diagnostic") != NULL) {
            p->diag->time_print = double_reader(line, len);
        } else if(strstr(line, "time_stop") != NULL) {
            p->time->time_stop = double_reader(line, len);
        } else if(strstr(line, "sigma_boltzman") != NULL) {
            p->constants->sigma_boltzman = double_reader(line, len);
        } else if(strstr(line, "source") != NULL) {
            p->constants->source = double_reader(line, len);
        } else if(strstr(line, "c_light") != NULL) {
            p->constants->c_light = double_reader(line, len);
        } else if(strstr(line, "sigma_factor") != NULL) {
            p->constants->sigma_factor = double_reader(line, len);
        } else if(strstr(line, "T0") != NULL) {
            p->constants->T0 = double_reader(line, len);
        } else if(strstr(line, "TH") != NULL) {
            p->constants->TH = double_reader(line, len);
        } else if(strstr(line, "bc_type") != NULL) {
            printf("??");
        } else if(strstr(line, "num_mats") != NULL) {
            p->mats->num_mats = double_reader(line, len);
        } else if(strstr(line, "mat_type") != NULL) {
            p->mats->mat_type = int_reader(line, len);
        } else if(strstr(line, "xmax") != NULL) {
            xmax = int_reader(line, len);
        }else if(strstr(line, "ymax") != NULL) {
            ymax = int_reader(line, len);
        } else if (strstr(line, "sig_fac") != NULL) {
            p->constants->sigma_factor = double_reader(line, len);
        }
    }
    p->constants->TH = p->constants->TH * HEV;
    p->coor->nxp = p->vol->nx + 1;
    p->coor->nyp = p->vol->ny + 1;
    p->constants->dx = xmax / (p->coor->nxp + 2) ;
    p->constants->dy = ymax / (p->coor->nyp + 2) ;

    fclose(fp);
}

void init_materials_datafile(Material *m, int mat_number) {
    FILE *fp;
    char * line = NULL;
   if (mat_number == SUOLSON) {
        fp = fopen("Su_olson", "r");
   }
   else if (mat_number == OLSON) {
       fp = fopen("olson", "r");
   }
    size_t len = 0;
    ssize_t read;
    while ((read = getline(&line, &len, fp)) != -1) {
        // set up the constants of the problem
        if (strstr(line, "alpha") != NULL) {
            m->alpha = double_reader(line, len);
        } else if(strstr(line, "lambda") != NULL) {
            m->lambda = double_reader(line, len);
        } else if(strstr(line, "beta") != NULL) {
            m->beta = double_reader(line, len);
        } else if(strstr(line, "mu") != NULL) {
            m->mu = double_reader(line, len);
        } else if(strstr(line, "g") != NULL) {
            m->g = double_reader(line, len);
        } else if(strstr(line, "f") != NULL) {
            m->f = double_reader(line, len);
        } else if(strstr(line, "i_start") != NULL) {
            m->i_start = int_reader(line, len);
        } else if(strstr(line, "i_end") != NULL) {
            m->i_end = int_reader(line, len) + 2;
        } else if(strstr(line, "j_start") != NULL) {
            m->j_start = int_reader(line, len);
        } else if(strstr(line, "j_end") != NULL) {
            m->j_end = int_reader(line, len) + 2;
            
        } else if(strstr(line, "rho") != NULL) {
            m->init_rho = double_reader(line, len);
        }
    }
    //normalize g and f
//f = 5.485E10
//תנאי שפה ב-hev
    m->g = m->g / pow(HEV, m->alpha);
    m->f = m->f / pow(HEV, m->beta);
    
}


/**
 * @brief initializes all of the structures
 * First it goes to another function that initialize a python interpeter
 * and executes the python function that will give us the sizes.
 * Second, we increment the number of K,L by 2, this is the imaginary cells.
 * Third, we allocate the memory.
 * Fourth initialize values.
 */
void init(Problem *p, char *datafaile) {
    int i, j;
    int nxp,nyp,nx,ny;
    //init_python(p);
    //first we malloc all structures

    init_datafile(p, datafaile);
    p->mats->mat = (Material *) malloc(sizeof(Material) * p->mats->num_mats);
    for (i = 0; i < p->mats->num_mats; i++) {
        init_materials_datafile(&p->mats->mat[i], p->mats->mat_type);
    }
    if (p->mats->num_mats == 1) {
        if (VERBOSE) 
            printf("One Material\n");
        p->mats->mat[0].i_start = 0;
        p->mats->mat[0].j_start = 0;
        p->mats->mat[0].i_end = p->vol->nx + 2; 
        p->mats->mat[0].j_end = p->vol->ny + 2; 
        
    }
    
    printf("Done reading Materials\n");
    
    // WE HAVE IMAGINARY CELLS. SO WE NEED TO ADD FOR THE VERTEX QUANT..
    // + 2 (right and left edge)
    // AND TO CELL QUANTITY + 2 ASWELL (EACH)
    nxp = p->coor->nxp += 2;
    nyp = p->coor->nyp += 2;
    // VOLUME NEEDS TO BE LAST!!!!!!!!!!!
    nx = p->diff_coeff->nx = p->energy->nx = p->temp->nx = p->rho->nx = p->heat_cap->nx = p->opacity->nx = p->vol->nx += 2;
    ny = p->diff_coeff->ny = p->energy->ny = p->temp->ny = p->rho->ny = p->heat_cap->ny = p->opacity->ny = p->vol->ny += 2;
    printf("%d %d\n",nx,ny);
    p->constants->a_rad = 4.0 * p->constants->sigma_boltzman / p->constants->c_light;
    //malloc
    p->coor->R = malloc_2d(nxp, nyp);
    p->coor->Z = malloc_2d(nxp, nyp);
    p->coor->X = malloc_2d(nxp, nyp);
    p->coor->Y = malloc_2d(nxp, nyp);
    
    p->energy->current = malloc_2d(nx, ny);
    p->energy->prev    = malloc_2d(nx, ny);
    p->temp->current   = malloc_2d(nx, ny);
    p->temp->prev      = malloc_2d(nx, ny);

    p->vol->values = malloc_2d(nx, ny);
    p->rho->values = malloc_2d(nx, ny);
    
    p->diff_coeff->values = malloc_2d(nx, ny);
    p->opacity->values    = malloc_2d(nx, ny);
    p->heat_cap->values   = malloc_2d(nx, ny);

    //time
    p->time->cycle = 0;
    p->time->time_passed = p->time->t0;
    //init
    // init_mesh_Kershaw1(p->coor->nxp, p->coor->nyp,p->coor->X, p->coor->Y, p->constants->dx);
    
    // for 1d problem on y-axis
    if (nx == 3) {
        for ( i = 0; i < p->coor->nxp; i++) {
            for ( j = 0; j < p->coor->nyp; j++) {
                p->coor->X[i][j] = p->coor->R[i][j] = 1;
                p->coor->Y[i][j] = p->coor->Z[i][j] = j * p->constants->dy;
            }
        }
        for ( i = 0; i < p->coor->nxp; i++) {
            for ( j = 0; j < p->coor->nyp; j++) {
                if (i == 0) {
                    // p->coor->X[i][j] = p->coor->R[i][j] = -1;
                    p->coor->X[i][j] = -1;
                }else if(i == 1) {
                    // p->coor->X[i][j] = p->coor->R[i][j] = 0;
                    p->coor->X[i][j] = 0;
                }else if(i == 2) {
                    p->coor->X[i][j] = 1;
                    // p->coor->X[i][j] = p->coor->R[i][j]= 1;
                }else if(i == 3) {
                    p->coor->X[i][j] = 2;
                    // p->coor->X[i][j] = p->coor->R[i][j]= 2;
                }
            }
        }
    }
    //for 1d problem on x-axis
    if (ny == 3) {
        for ( i = 0; i < p->coor->nxp; i++) {
            for ( j = 0; j < p->coor->nyp; j++) {
                // p->coor->X[i][j] = p->coor->R[i][j] = 1;
                p->coor->X[i][j] = i * p->constants->dx;
                p->coor->R[i][j] = 1;
                p->coor->Z[i][j] = 0;
                
            }
        }
        
        for ( i = 0; i < p->coor->nxp; i++) {
            for ( j = 0; j < p->coor->nyp; j++) {
                if (j == 0) {
                    p->coor->Y[i][j] = -1;
                }else if(j == 1) {
                    p->coor->Y[i][j] = 0;
                }else if(j == 2) {
                    p->coor->Y[i][j] =  1;
                }else if(j == 3) {
                    p->coor->Y[i][j] = 2;
                }
            }
        }
            // init_mesh_Kershaw1(p->coor->nxp, p->coor->nyp,p->coor->X, p->coor->Y, p->constants->dx);

    }

    //BC for X
    
    for ( i = 0 ; i < nx; i++) {
        for (j = 0 ; j < ny; j++) {
            p->energy->prev[i][j] = p->energy->current[i][j] = pow(65255,4) * p->constants->a_rad;//p->constants->a_rad * pow(62525.5, 4);
            p->temp->prev[i][j] = p->temp->current[i][j] = 65255;//*//10E-5 * HEV;//pow(p->constants->T0, 4) * p->constants->a_rad;
        }
    }
    if (ny == 3) {
        for (j = 0; j < ny; j++){
            p->energy->prev[0][j]    = p->constants->a_rad * pow(p->constants->TH, 4);
            p->energy->current[0][j] = p->constants->a_rad * pow(p->constants->TH, 4);
            p->temp->prev[0][j]      = p->constants->TH; //pow(p->constants->TH,4);// * p->constants->a_rad;
            p->temp->current[0][j]   = p->constants->TH; //pow(p->constants->TH,4);// * p->constants->a_rad;
        }
    }
    if (nx == 3) {
        for (i = 0; i < nx; i++){
            p->energy->prev[i][0]    = p->constants->a_rad * pow(p->constants->TH, 4);
            p->energy->current[i][0] = p->constants->a_rad * pow(p->constants->TH, 4);
            p->temp->prev[i][0]      = p->constants->TH; //pow(p->constants->TH,4);// * p->constants->a_rad;
            p->temp->current[i][0]   = p->constants->TH; //pow(p->constants->TH,4);// * p->constants->a_rad;
        }
    }
    
    mesh_square_volume(p->vol->values, p->coor->X,p->coor->Y, nx, ny);
    init_density(p->mats, p->rho);
    diagnostics_initial(p);
    printf("Done init\n");
}

/**
 * @brief Initailizes the mesh so that it will have constant delta R and delta Z
 * 
 * **/
void init_mesh_Kershaw1(int nxp, int nyp, double **R, double **Z, double dr) {
    int i = 0, j = 0;
    double dz = 0.1;
    for ( i = 0; i < nxp; i++) {
        for (j = 0 ; j < nyp; j++) {
            R[i][j] =  i * dr;

        }
    }
    
    for ( i = 0 ; i < nxp; i++) {
        for (j = 0 ; j < nyp; j++) {
            if (j == 0) {
                Z[i][j] = -1;
            }else if(j == 1) {
                Z[i][j] = 0;
            }else if(j == 2) {
                Z[i][j] = 1;
            }else if(j == 3) {
                Z[i][j] = 2;
            }
            
        }
    }
    
    // for ( i = 1 ; i < nxp - 1; i++) {
    //     for (j = 1 ; j < nyp - 1; j++) {
    //         Z[i][j] = 0;
    //     }
    // }
}

/*
* @brief Initalizes the density of the problem.
*
**/
void init_density(Materials *mats, Data *density) {
    int i_start, i_end, j_start, j_end, k;
    int i,j ;
    double init_rho;
    double **rho = density->values;
    for (k = 0; k < mats->num_mats; k++) {
        i_start = mats->mat[k].i_start;
        i_end = mats->mat[k].i_end;
        j_start = mats->mat[k].j_start;
        j_end = mats->mat[k].j_end;
        init_rho = mats->mat[k].init_rho;
        for (i = i_start; i < i_end; i++) {
            for (j = j_start; j < j_end; j++) {
                rho[i][j] = init_rho;
            }
        }
    }
}

void clean_prog(Problem *p) {
    int nxp = p->coor->nxp;
    int nyp = p->coor->nyp;
    int nx = p->energy->nx;
    int ny = p->energy->ny;
    free_2d(p->coor->R,nxp);
    free_2d(p->coor->Z,nxp);
    free_2d(p->energy->current,nx);
    free_2d(p->energy->prev,nx);
    free_2d(p->vol->values,nx);
    free(p->coor);
    free(p->diff_coeff);
    free(p->vol);
    free(p->rho);
    free(p->opacity);
    free(p->heat_cap);
    free(p->time);
    free(p->energy);
    //free_3d(A,nxp,nyp);
}
