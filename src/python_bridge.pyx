import module
from module import input
cimport cython

cdef public void setup_params():
    print "hello from cython"
    #input.print_py1()
    return;
