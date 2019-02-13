#include <Python.h>
#include "structmember.h"

typedef struct {
    PyObject_HEAD
    int k_max;
    int l_max;
    int kc_max;
    int lc_max;
    double dt;
    double t0;
    double time_diagnostic;
    double time_stop;
} Datafile;

static void
Datafile_dealloc(Datafile *self)
{
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
Datafile_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Datafile *self;
    self = (Datafile *) type->tp_alloc(type, 0);
    if (self != NULL) {
        self->k_max  = 0;
        self->l_max  = 0;
        self->kc_max = 0;
        self->lc_max = 0;
        self->time_stop = 0;
        self->dt = 0;
        self->t0 = 0;
        self->time_diagnostic = 0.0;
    }
    return (PyObject *) self;
}

static int
Datafile_init(Datafile *self, PyObject *args, PyObject *kwds) {
   static char *kwlist[] = {"first", "last", "number", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|iiiidddd", kwlist,
                                      &self->k_max,
                                      &self->l_max,
                                      &self->kc_max,
                                      &self->lc_max,
                                      &self->time_stop,
                                      &self->dt,
                                      &self->t0,
				                      &self->time_diagnostic
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
     {"dt", T_DOUBLE, offsetof(Datafile, dt), 0,
     "dt"},
     {"time_stop", T_DOUBLE, offsetof(Datafile, time_stop), 0,
     "time_stop"},
     {"t0", T_DOUBLE, offsetof(Datafile, t0), 0,
     "t0"},
     {"time_diagnostic", T_DOUBLE, offsetof(Datafile, time_diagnostic), 0,
     "time_diagnostic"},
    {NULL}  /* Sentinel */
};

static PyObject *
Datafile_name(Datafile *self, PyObject *Py_UNUSED(ignored)) {
    return NULL;
}

static PyMethodDef Datafile_methods[] = {
    {"name", (PyCFunction) Datafile_name, METH_NOARGS,
     "Return the name, combining the first and last name"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject DatafileType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "df.Datafile",
    .tp_doc = "Datafile objects",
    .tp_basicsize = sizeof(Datafile),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = Datafile_new,
    .tp_init = (initproc) Datafile_init,
    .tp_dealloc = (destructor) Datafile_dealloc,
    .tp_members = Datafile_members,
    .tp_methods = Datafile_methods,
};

static PyModuleDef datafilemodule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "df",
    .m_doc = "Example module that creates an extension type.",
    .m_size = -1,
};

PyMODINIT_FUNC
PyInit_df(void)
{
    PyObject *m;
    if (PyType_Ready(&DatafileType) < 0)
        return NULL;

    m = PyModule_Create(&datafilemodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&DatafileType);
    PyModule_AddObject(m, "Datafile", (PyObject *) &DatafileType);
    return m;
}