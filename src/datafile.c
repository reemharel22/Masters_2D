#include <Python.h>
#include "python2.7/structmember.h"

typedef struct {
    PyObject_HEAD
    int k_max;
    int l_max;
    int kc_max;
    int lc_max;
    int bc_type;
    int num_mats;
    int i_start;
    int i_end;
    int j_start;
    int j_end;
    double dt;
    double t0;
    double time_diagnostic;
    double time_stop;
    double a_rad;
    double sigma_boltzman;
    double c;
    double g;
    double f;
    double dt_max;
    double alpha;
    double beta;
    double lambda1;
    double rho;
    double mu;
    double T0;
    double dt_factor;
} Datafile;

static void
Datafile_dealloc(Datafile* self)
{
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
Datafile_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Datafile *self;

    self = (Datafile *)type->tp_alloc(type, 0);
    self->k_max  = 0;
    self->l_max  = 0;
    self->kc_max = 0;
    self->lc_max = 0;
    self->num_mats = 0;
    self->i_start = 0;
    self->i_end = 0;
    self->j_start = 0;
    self->j_end = 0;
    self->bc_type = 0;
    self->time_stop = 0;
    self->dt = 0;
    self->t0 = 0;
    self->time_diagnostic = 0.0;
    self->a_rad = 0.0;
    self->sigma_boltzman = 0.0;
    self->c = 0.0;
    self->g = 0.0;
    self->f = 0.0;
    self->alpha = 0.0;
    self->beta = 0.0;
    self->lambda1 = 0.0;
    self->mu = 0.0;
    self->rho = 0.0;
    self->T0 = 0.0;
    self->dt_factor = 0.0;
    self->dt_max = 0.0;

    return (PyObject *)self;
}

static int
Datafile_init(Datafile *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"first", "last", "number", NULL};
    
    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|iiiiiiiiiidddddddddddddddd", kwlist,
                                        &self->k_max,
                                        &self->l_max,
                                        &self->kc_max,
                                        &self->lc_max,
                                        &self->num_mats,
                                        &self->i_start,
                                        &self->i_end,
                                        &self->j_start,
                                        &self->j_end,
                                        &self->bc_type,
                                        &self->time_stop,
                                        &self->dt,
                                        &self->t0,
                                        &self->time_diagnostic,
                                        &self->a_rad,
                                        &self->sigma_boltzman,
                                        &self->c,
                                        &self->g,
                                        &self->f,
                                        &self->alpha,
                                        &self->beta,
                                        &self->lambda1,
                                        &self->mu,
                                        &self->rho,
                                        &self->T0,
                                        &self->dt_factor,
                                        &self->dt_max
                                      ))
        return -1;


    return 0;
}


static PyMemberDef Datafile_members[] = {
    {"k_max", T_INT, offsetof(Datafile, k_max), 0,
     "k_max"},
    {"l_max", T_INT, offsetof(Datafile, l_max), 0,
     "l_max"},
     {"kc_max", T_INT, offsetof(Datafile, kc_max), 0,
     "xx"},
    {"lc_max", T_INT, offsetof(Datafile, lc_max), 0,
     "lc_max"},
      {"bc_type", T_INT, offsetof(Datafile, bc_type), 0,
     "bc_type"},
     {"num_mats", T_INT, offsetof(Datafile, num_mats), 0,
     "num_mats"},
      {"i_start", T_INT, offsetof(Datafile, i_start), 0,
     "i_start"},
      {"i_end", T_INT, offsetof(Datafile, i_end), 0,
     "i_end"},
      {"j_start", T_INT, offsetof(Datafile, j_start), 0,
     "j_start"},
      {"j_end", T_INT, offsetof(Datafile, j_end), 0,
     "j_end"},
     {"dt", T_DOUBLE, offsetof(Datafile, dt), 0,
     "dt"},
     {"time_stop", T_DOUBLE, offsetof(Datafile, time_stop), 0,
     "time_stop"},
     {"t0", T_DOUBLE, offsetof(Datafile, t0), 0,
     "t0"},
     {"time_diagnostic", T_DOUBLE, offsetof(Datafile, time_diagnostic), 0,
     "time_diagnostic"},
     {"a_rad", T_DOUBLE, offsetof(Datafile, a_rad), 0,
     "a_rad"},
     {"sigma_boltzman", T_DOUBLE, offsetof(Datafile, sigma_boltzman), 0,
     "sigma_boltzman"},
     {"c", T_DOUBLE, offsetof(Datafile, c), 0,
     "c"},
     {"g", T_DOUBLE, offsetof(Datafile, g), 0,
     "g"},
     {"f", T_DOUBLE, offsetof(Datafile, f), 0,
     "f"},
     {"alpha", T_DOUBLE, offsetof(Datafile, alpha), 0,
     "alpha"},
     {"beta", T_DOUBLE, offsetof(Datafile, beta), 0,
     "beta"},
     {"lambda1", T_DOUBLE, offsetof(Datafile, lambda1), 0,
     "lambda1"},
     {"mu", T_DOUBLE, offsetof(Datafile, mu), 0,
     "mu"},
      {"rho", T_DOUBLE, offsetof(Datafile, rho), 0,
     "rho"},
     {"T0", T_DOUBLE, offsetof(Datafile, T0), 0,
     "T0"},
     {"dt_factor", T_DOUBLE, offsetof(Datafile, dt_factor), 0,
     "dt_factor"},
     {"dt_max", T_DOUBLE, offsetof(Datafile, dt_max), 0,
     "dt_max"},
    {NULL}  /* Sentinel */
};

static PyObject *
Datafile_name(Datafile* self)
{

    return NULL;
}

static PyMethodDef Datafile_methods[] = {
    {"name", (PyCFunction)Datafile_name, METH_NOARGS,
     "Return the name, combining the first and last name"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject DatafileType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "datafile.Datafile",             /* tp_name */
    sizeof(Datafile),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Datafile_dealloc, /* tp_dealloc */
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
    "Datafile objects",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Datafile_methods,             /* tp_methods */
    Datafile_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Datafile_init,      /* tp_init */
    0,                         /* tp_alloc */
    Datafile_new,                 /* tp_new */
};

static PyMethodDef module_methods1[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initdf(void)
{
    PyObject* m;

    if (PyType_Ready(&DatafileType) < 0)
        return;

    m = Py_InitModule3("df", module_methods1,
                       "Example module that creates an extension type.");

    if (m == NULL)
        return;

    Py_INCREF(&DatafileType);
    PyModule_AddObject(m, "Datafile", (PyObject *)&DatafileType);
}
