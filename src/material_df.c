#include <Python.h>
#include "python2.7/structmember.h"

typedef struct {
    PyObject_HEAD
    int num_mats;
    int i_start;
    int i_end;
    int j_start;
    int j_end;
    double g;
    double f;
    double alpha;
    double beta;
    double lambda1;
    double mu;
} Materialdf;

static void
Materialdf_dealloc(Materialdf* self)
{
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
Materialdf_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Materialdf *self;

    self = (Materialdf *)type->tp_alloc(type, 0);
    self->i_start = 0;
    self->i_end = 0;
    self->j_start = 0;
    self->j_end = 0;
    self->g = 0.0;
    self->f = 0.0;
    self->alpha = 0.0;
    self->beta = 0.0;
    self->lambda1 = 0.0;
    self->mu = 0.0;

    return (PyObject *)self;
}

static int
Materialdf_init(Materialdf *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"first", "last", "number", NULL};
    
    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|iiiidddddd", kwlist,
                                        &self->i_start,
                                        &self->i_end,
                                        &self->j_start,
                                        &self->j_end,
                                        &self->g,
                                        &self->f,
                                        &self->alpha,
                                        &self->beta,
                                        &self->lambda1,
                                        &self->mu
                                      ))
        return -1;


    return 0;
}


static PyMemberDef Materialdf_members[] = {
    {"k_max", T_INT, offsetof(Materialdf, k_max), 0,
     "k_max"},
    {"l_max", T_INT, offsetof(Materialdf, l_max), 0,
     "l_max"},
     {"kc_max", T_INT, offsetof(Materialdf, kc_max), 0,
     "xx"},
    {"lc_max", T_INT, offsetof(Materialdf, lc_max), 0,
     "lc_max"},
      {"bc_type", T_INT, offsetof(Materialdf, bc_type), 0,
     "bc_type"},
     {"num_mats", T_INT, offsetof(Materialdf, num_mats), 0,
     "num_mats"},
      {"i_start", T_INT, offsetof(Materialdf, i_start), 0,
     "i_start"},
      {"i_end", T_INT, offsetof(Materialdf, i_end), 0,
     "i_end"},
      {"j_start", T_INT, offsetof(Materialdf, j_start), 0,
     "j_start"},
      {"j_end", T_INT, offsetof(Materialdf, j_end), 0,
     "j_end"},
     {"dt", T_DOUBLE, offsetof(Materialdf, dt), 0,
     "dt"},
     {"time_stop", T_DOUBLE, offsetof(Materialdf, time_stop), 0,
     "time_stop"},
     {"t0", T_DOUBLE, offsetof(Materialdf, t0), 0,
     "t0"},
     {"time_diagnostic", T_DOUBLE, offsetof(Materialdf, time_diagnostic), 0,
     "time_diagnostic"},
     {"a_rad", T_DOUBLE, offsetof(Materialdf, a_rad), 0,
     "a_rad"},
     {"sigma_boltzman", T_DOUBLE, offsetof(Materialdf, sigma_boltzman), 0,
     "sigma_boltzman"},
     {"c", T_DOUBLE, offsetof(Materialdf, c), 0,
     "c"},
     {"g", T_DOUBLE, offsetof(Materialdf, g), 0,
     "g"},
     {"f", T_DOUBLE, offsetof(Materialdf, f), 0,
     "f"},
     {"alpha", T_DOUBLE, offsetof(Materialdf, alpha), 0,
     "alpha"},
     {"beta", T_DOUBLE, offsetof(Materialdf, beta), 0,
     "beta"},
     {"lambda1", T_DOUBLE, offsetof(Materialdf, lambda1), 0,
     "lambda1"},
     {"mu", T_DOUBLE, offsetof(Materialdf, mu), 0,
     "mu"},
     {"T0", T_DOUBLE, offsetof(Materialdf, T0), 0,
     "T0"},
     {"dt_factor", T_DOUBLE, offsetof(Materialdf, dt_factor), 0,
     "dt_factor"},
     {"dt_max", T_DOUBLE, offsetof(Materialdf, dt_max), 0,
     "dt_max"},
    {NULL}  /* Sentinel */
};

static PyObject *
Materialdf_name(Materialdf* self)
{

    return NULL;
}

static PyMethodDef Materialdf_methods[] = {
    {"name", (PyCFunction)Materialdf_name, METH_NOARGS,
     "Return the name, combining the first and last name"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject MaterialdfType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Materialdf.Materialdf",             /* tp_name */
    sizeof(Materialdf),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Materialdf_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_compare */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Materialdf objects",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Materialdf_methods,             /* tp_methods */
    Materialdf_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Materialdf_init,      /* tp_init */
    0,                         /* tp_alloc */
    Materialdf_new,                 /* tp_new */
};

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initdf(void)
{
    PyObject* m;

    if (PyType_Ready(&MaterialdfType) < 0)
        return;

    m = Py_InitModule3("df", module_methods,
                       "Example module that creates an extension type.");

    if (m == NULL)
        return;

    Py_INCREF(&MaterialdfType);
    PyModule_AddObject(m, "Materialdf", (PyObject *)&MaterialdfType);
}
