#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "initialize.h"
//#include "python2.7/Python.h"
#include "diagnostics.h"
//#include "datafile.c"
/**
 * @file initialize.c
 */

/**
 * @brief initializes the structs.
 */
void init_python(Problem*p) {
   /* PyObject *pName, *pModule, *pFunc;
    char *argv[3] = {"call","input","test"};
    Datafile *n = 0;

    Py_Initialize();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    pName = PyString_FromString(argv[1]);
    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, argv[2]);
        /* pFunc is a new reference */
    
      /*  if (pFunc && PyCallable_Check(pFunc)) {
           /* pArgs = PyTuple_New(argc - 3);
            for (i = 0; i < argc - 3; ++i) {
                pValue = PyInt_FromLong(atoi(argv[i + 3]));
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    return 1;
                }
                PyTuple_SetItem(pArgs, i, pValue);
            }*/

      /*      n =(Datafile*) PyObject_CallObject(pFunc, 0);
            //pValue = PyObject_CallObject(pFunc, pArgs);
            //PyArg_ParseTupleAndKeywords(Pva)
            if (n != NULL) {
                p->coor.nxp = n->nxp;
                p->coor.nyp = n->nyp;
                p->vol.nx =  n->nx;
                p->vol.ny =  n->ny;
                p->time.t0 = n->t0;
                p->time.dt = n->dt;
                p->time.time_stop = n->time_stop;
                p->diag.time_print = n->time_diagnostic;
                p->constants.a_rad = n->a_rad;
                p->constants.sigma_boltzman = n->sigma_boltzman;
                p->constants.c_light = n->c;
                p->constants.T0 = n->T0;
                p->time.dt_factor = n->dt_factor;
                p->boundary_type = n->bc_type;
                p->mats.num_mats = n->num_mats;
                Datafile_dealloc(n);
                Py_DECREF(n);
               // Py_DECREF(pArgs);

            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                exit(1);
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", argv[2]);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", argv[1]);
        Py_Finalize();  
        exit(1);
    }
    Py_Finalize();*/
}

/*
* @brief Initializes the materials thro python
*/
void init_materials_python(Material *m, int mat_number) {
  /*  PyObject *pName, *pModule, *pFunc;
    char* argv[3] = {"call", "input", "material1"}; 
    if (mat_number == 1) {
        argv[2] = "material2";    
    }
    Datafile *n = 0;
    Py_Initialize();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    pName = PyString_FromString(argv[1]);
    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, argv[2]);
        /* pFunc is a new reference */
    
      /*  if (pFunc && PyCallable_Check(pFunc)) {
           /* pArgs = PyTuple_New(argc - 3);
            for (i = 0; i < argc - 3; ++i) {
                pValue = PyInt_FromLong(atoi(argv[i + 3]));
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    return 1;
                }
                PyTuple_SetItem(pArgs, i, pValue);
            }*/

      /*      n = (Datafile*) PyObject_CallObject(pFunc, 0);
            //pValue = PyObject_CallObject(pFunc, pArgs);
            //PyArg_ParseTupleAndKeywords(Pva)
            if (n != NULL) {
                m->alpha = n->alpha; // !< ALpha, related to kappa rossland
                m->lambda = n->lambda1; // !< Lambda related to cv
                m->beta = n->beta; // !< BEta related to rossaland
                m->mu = n->mu; // !< Mu related to rossland
                m->g = n->g;  //!< g Related to rossland
                m->f = n->f; // !< f related to rossland
                m->i_start = n->i_start;
                m->i_end = n->i_end;
                m->j_start = n->j_start;
                m->j_end = n->j_end;
                m->init_rho = n->rho;
                Datafile_dealloc(n);
                Py_DECREF(n);
               // Py_DECREF(pArgs);

            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                exit(1);
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", argv[2]);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", argv[1]);
        Py_Finalize();  
        exit(1);
    }
    
    Py_Finalize();*/
}

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
        } else if(strstr(line, "c") != NULL) {
            p->constants->c_light = double_reader(line, len);
        } else if(strstr(line, "sigma_factor") != NULL) {
            p->constants->sigma_factor = double_reader(line, len);
        } else if(strstr(line, "T0") != NULL) {
            p->constants->T0 = double_reader(line, len);
        } else if(strstr(line, "bc_type") != NULL) {
            printf("??");
        } else if(strstr(line, "num_mats") != NULL) {
            p->mats->num_mats = int_reader(line, len);
        } else if(strstr(line, "mat_type") != NULL) {
            p->mats->mat_type = int_reader(line, len);
        } else if (strstr(line, "sig_fac") != NULL) {
                    printf("A!Q@W#E$RTYG\n");

            p->constants->sigma_factor = double_reader(line, len);
        }
    }
    
    fclose(fp);
}

void init_materials_datafile(Material *m, int mat_number) {
    FILE *fp;
    char * line = NULL;
  //  if (mat_number == 1) {
        fp = fopen("Silicon_Dioxide", "r");//SILICON DIOXIDE
   // }
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
    m->g = m->g / pow(1160500.0, m->alpha);
    m->f = m->f / pow(1160500.0, m->beta);
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
    printf("Done reading Materials\n");
    
    // WE HAVE IMAGINARY CELLS. SO WE NEED TO ADD FOR THE VERTEX QUANT..
    // + 2 (right and left edge)
    // AND TO CELL QUANTITY + 2 ASWELL (EACH)
    nxp = p->coor->nxp += 2;
    nyp = p->coor->nyp += 2;
    nx = p->diff_coeff->nx = p->energy->nx = p->temp->nx = p->vol->nx += 2;
    ny = p->diff_coeff->ny = p->energy->ny = p->temp->ny = p->vol->ny += 2;

    p->constants->a_rad = 4.0 * p->constants->sigma_boltzman / p->constants->c_light;
    //malloc
    p->coor->R = malloc_2d(nxp, nyp );
    p->coor->Z = malloc_2d(nxp, nyp );
    
    
    p->energy->current = malloc_2d(nx, ny );
    p->energy->prev    = malloc_2d(nx, ny );
    p->temp->current   = malloc_2d(nx, ny );
    p->temp->prev      = malloc_2d(nx, ny );

    p->vol->values = malloc_2d(nx, ny);
    p->rho->values = malloc_2d(nx, ny);
    
    p->diff_coeff->values = malloc_2d(nx, ny);
    p->opacity->values    = malloc_2d(nx, ny);
    p->heat_cap->values   = malloc_2d(nx, ny);

    //time
    p->time->cycle = 0;
    p->time->time_passed = p->time->t0;
    
    //init
    init_mesh_Kershaw1(p->coor->nxp, p->coor->nyp,p->coor->R, p->coor->Z);
    for ( i = 0 ; i < nx; i++) {
        for (j = 0 ; j < ny; j++) {
            p->energy->prev[i][j] = p->energy->current[i][j] = p->temp->prev[i][j] = p->temp->current[i][j] = pow(p->constants->T0, 4) * p->constants->a_rad;
        }
    }
    
    mesh_square_volume(p->vol->values, p->coor->R,p->coor->Z, nx, ny);
    init_density(p->mats, p->rho);
    diagnostics_initial(p);
}

/**
 * @brief Initailizes the mesh so that it will have constant delta R and delta Z
 * 
 * **/
void init_mesh_Kershaw1(int nxp, int nyp, double **R, double **Z) {
    int i = 0, j = 0;
    double dr = 0.01;
    double dz = 0.01;
    for ( i = 0; i < nxp; i++) {
        for (j = 0 ; j < nyp; j++) {
            R[i][j] = j * dr;
        }
    }
    
    for ( i = 0 ; i < nxp; i++) {
        for (j = 0 ; j < nyp; j++) {
            Z[i][j] = i * dz;
        }
    }
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
