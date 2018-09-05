#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "initialize.h"
#include "python3.6m/Python.h"
#include "python_bridge.h"

/**
 * @file initialize.c
 */

/**
 * @brief initializes the structs.
 */
void init1(Problem*p) {
    PyObject *pName, *pModule, *pFunc;
    PyObject *pArgs, *pValue;
    PyObject *ob1,*ob2;
    int i;
    int argc = 5;
    double a = 2;
    
    char*argv[5] = {"call","input","multiply","4","5"};
    
    Py_Initialize();
    
    //PySys_SetArgv(argc,argv);
    setup_params();
    printf("hello\n");
    Py_Finalize();
    return;
    /*pName = PyString_FromString(argv[1]);
    /* Error checking of pName left out 

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, argv[2]);
        /* pFunc is a new reference 

        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(argc - 3);
            for (i = 0; i < argc - 3; ++i) {
                pValue = PyInt_FromLong(atoi(argv[i + 3]));
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    return 1;
                }
                PyTuple_SetItem(pArgs, i, pValue);
            }
            pValue = PyEval_CallObject(pFunc, pArgs);
            //pValue = PyObject_CallObject(pFunc, pArgs);
            //PyArg_ParseTupleAndKeywords(Pva)
            if (pValue != NULL) {
                PyArg_ParseTuple(pValue,"[i]",ob1,ob2);
                printf("Result of call: %ld\t%ld\n", PyInt_AsLong(pValue),PyInt_AsLong(ob2));
                printf("hello\n");
                exit(1);
                Py_DECREF(pValue);
                Py_DECREF(pArgs);

            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
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
        return 1;
    }
    Py_Finalize();*/
}
void init(Problem*p) {
    int i,j;
    double K_max,L_max, KC_max,LC_max;
    double **R,**Z,**E,**E_o;
    //coordinates 
    K_max = 78;
    L_max = 76;
    KC_max = K_max - 1;
    LC_max = L_max - 1;
    // WE HAVE IMAGINARY CELLS. SO WE NEED TO ADD FOR THE VERTEX QUANT..
    // + 2 (right and left edge)
    // AND TO CELL QUANTITY + 2 ASWELL (EACH)
    p->coor.K_max = K_max;
    p->coor.L_max = L_max;
    p->coor.R = malloc_2d(K_max, L_max );
    p->coor.Z = malloc_2d(K_max, L_max );
    
    //energy 
    p->eng.KC_max = p->vol.KC_max =  KC_max;
    p->eng.LC_max = p->vol.LC_max =  LC_max;
    p->eng.E_current = malloc_2d(KC_max, LC_max );
    p->eng.E_old     = malloc_2d(KC_max, LC_max );

    //volume
    p->vol.volume = malloc_2d(KC_max, LC_max);

    //time
    p->time.cycle = 0;
    p->time.t0 = 0;
    p->time.dt = 0.01;
    p->time.time_passed = p->time.t0;
    p->time.time_stop = 0.01;

    //init
    init_mesh_Kershaw1(p->coor.K_max,p->coor.L_max,p->coor.R,p->coor.Z);

    for ( i = 0 ; i < KC_max; i++) {
        for (j = 0 ; j < LC_max; j++) {
             p->eng.E_current[i][j] = 0;
        }
    }
    mesh_square_volume(p->vol.volume, p->coor.R,p->coor.Z,K_max,L_max);
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
    int KC_max = p->eng.KC_max;
    int LC_max = p->eng.LC_max;
    free_2d(p->coor.R,K_max);
    free_2d(p->coor.Z,K_max);
    free_2d(p->eng.E_current,KC_max);
    free_2d(p->eng.E_old,KC_max);
    free_2d(p->vol.volume,KC_max);
    //free_3d(A,K_max,L_max);
}
