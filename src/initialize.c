#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "initialize.h"
#include "python2.7/Python.h"
#include "datafile.c"
/**
 * @file initialize.c
 */

/**
 * @brief initializes the structs.
 */
void init_python(Problem*p) {
    PyObject *pName, *pModule, *pFunc;
    char *argv[3] = {"call","input","test"};
    Datafile *n = 0;

    /* Error checking of pName left out */
    Py_Initialize();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    pName = PyString_FromString(argv[1]);
    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, argv[2]);
        /* pFunc is a new reference */
    
        if (pFunc && PyCallable_Check(pFunc)) {
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

            n =(Datafile*) PyObject_CallObject(pFunc, 0);
            //pValue = PyObject_CallObject(pFunc, pArgs);
            //PyArg_ParseTupleAndKeywords(Pva)
            if (n != NULL) {
                p->coor.K_max = n->k_max;
                p->coor.L_max = n->l_max;
                p->vol.KC_max =  n->kc_max;
                p->vol.LC_max =  n->lc_max;
                p->time.t0 = n->t0;
                p->time.dt = n->dt;
                p->time.time_stop = n->time_stop;
                p->diag.time_print = n->time_diagnostic;
                p->constants.a_rad = n->a_rad;
                p->constants.sigma_boltzman = n->sigma_boltzman;
                p->constants.c_light = n->c;
                p->constants.g = n->g;
                p->constants.f = n->f; 
                p->constants.alpha = n->alpha;
                p->constants.beta = n->beta;
                p->constants.lambda = n->lambda1;
                p->constants.mu = n->mu; 
                p->constants.T0 = n->T0;
                p->boundary_type = n->bc_type;
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
    Py_Finalize();
}

/**
 * @brief initializes all of the structures
 * First it goes to another function that initialize a python interpeter
 * and executes the python function that will give us the sizes.
 * Second, we increment the number of K,L by 2, this is the imaginary cells.
 * Third, we allocate the memory.
 * Fourth initialize values.
 */
void init(Problem*p) {
    int i,j;
    int K_max,L_max,KC_max,LC_max;
    init_python(p);
    // WE HAVE IMAGINARY CELLS. SO WE NEED TO ADD FOR THE VERTEX QUANT..
    // + 2 (right and left edge)
    // AND TO CELL QUANTITY + 2 ASWELL (EACH)
    K_max = p->coor.K_max += 2;
    L_max = p->coor.L_max += 2;
    KC_max = p->diff_coeff.KC_max = p->energy.KC_max = p->vol.KC_max += 2;
    LC_max = p->diff_coeff.LC_max = p->energy.LC_max = p->vol.LC_max += 2;

    //malloc
    p->coor.R = malloc_2d(K_max, L_max );
    p->coor.Z = malloc_2d(K_max, L_max );

    p->energy.current = malloc_2d(KC_max, LC_max );
    p->energy.prev    = malloc_2d(KC_max, LC_max );
    
    p->vol.values = malloc_2d(KC_max, LC_max);
    p->rho.values = malloc_2d(KC_max, LC_max);
    p->diff_coeff.values = malloc_2d(KC_max, LC_max);

    //time
    p->time.cycle = 0;
    p->time.time_passed = p->time.t0;

    printf("Init done, with the following parameters:\n");
    printf("Time finish: %10e. Time start: %10e\n", p->time.time_passed, p->time.t0);
    printf("Diagnostic every: %10e time passed\n", p->diag.time_print);

    //init
    init_mesh_Kershaw1(p->coor.K_max,p->coor.L_max,p->coor.R,p->coor.Z);

    for ( i = 0 ; i < KC_max; i++) {
        for (j = 0 ; j < LC_max; j++) {
             p->energy.prev[i][j] = p->energy.current[i][j] = 0;
             p->temp.prev[i][j] = p->temp.current[i][j] = p->constants.T0;
        }
    }

    mesh_square_volume(p->vol.values, p->coor.R,p->coor.Z,KC_max,LC_max);

}

void init_mesh_Kershaw1(int K_max, int L_max, double **R, double **Z) {
    int i = 0, j = 0;
    for ( i = 0; i < K_max; i++) {
        for (j = 0 ; j < L_max; j++) {
            R[i][j] = (double) j / (L_max - 1);
        }
    }
    for ( i = 0 ; i < K_max; i++) {
        for (j = 0 ; j < L_max; j++) {
            Z[i][j] = (double) i / (K_max - 1);
        }
    }
}

void clean_prog(Problem *p) {
    int K_max = p->coor.K_max;
    int L_max = p->coor.L_max;
    int KC_max = p->energy.KC_max;
    int LC_max = p->energy.LC_max;
    free_2d(p->coor.R,K_max);
    free_2d(p->coor.Z,K_max);
    free_2d(p->energy.current,KC_max);
    free_2d(p->energy.prev,KC_max);
    free_2d(p->vol.values,KC_max);
    //free_3d(A,K_max,L_max);
}
