#include <Python.h>
#include "python2.7/structmember.h"

typedef struct {
    PyObject_HEAD
    int k_max;
    int l_max;
    int kc_max;
    int lc_max;
    double dt;
    double t0;
    double time_stop;
} Noddy;

static void
Noddy_dealloc(Noddy* self)
{
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
Noddy_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Noddy *self;

    self = (Noddy *)type->tp_alloc(type, 0);
    self->k_max  = 0;
    self->l_max  = 0;
    self->kc_max = 0;
    self->lc_max = 0;
    self->time_stop =0;
    self->dt = 0;
    self->t0 = 0;

    return (PyObject *)self;
}

static int
Noddy_init(Noddy *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"first", "last", "number", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|iiiiddd", kwlist,
                                      &self->k_max,
                                      &self->l_max,
                                      &self->kc_max,
                                      &self->lc_max,
                                      &self->time_stop,
                                      &self->dt,
                                      &self->t0
                                      ))
        return -1;


    return 0;
}


static PyMemberDef Noddy_members[] = {
    {"k_max", T_INT, offsetof(Noddy, k_max), 0,
     "k_max"},
    {"l_max", T_INT, offsetof(Noddy, l_max), 0,
     "l_max"},
     {"kc_max", T_INT, offsetof(Noddy, kc_max), 0,
     "xx"},
    {"lc_max", T_INT, offsetof(Noddy, lc_max), 0,
     "lc_max"},
     {"dt", T_DOUBLE, offsetof(Noddy, dt), 0,
     "dt"},
     {"time_stop", T_DOUBLE, offsetof(Noddy, time_stop), 0,
     "time_stop"},
     {"t0", T_DOUBLE, offsetof(Noddy, t0), 0,
     "t0"},
    {NULL}  /* Sentinel */
};

static PyObject *
Noddy_name(Noddy* self)
{

    return NULL;
}

static PyMethodDef Noddy_methods[] = {
    {"name", (PyCFunction)Noddy_name, METH_NOARGS,
     "Return the name, combining the first and last name"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject NoddyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "noddy.Noddy",             /* tp_name */
    sizeof(Noddy),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Noddy_dealloc, /* tp_dealloc */
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
    "Noddy objects",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Noddy_methods,             /* tp_methods */
    Noddy_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Noddy_init,      /* tp_init */
    0,                         /* tp_alloc */
    Noddy_new,                 /* tp_new */
};

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initnoddy2(void)
{
    PyObject* m;

    if (PyType_Ready(&NoddyType) < 0)
        return;

    m = Py_InitModule3("noddy2", module_methods,
                       "Example module that creates an extension type.");

    if (m == NULL)
        return;

    Py_INCREF(&NoddyType);
    PyModule_AddObject(m, "Noddy", (PyObject *)&NoddyType);
}