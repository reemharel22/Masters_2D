#!/bin/bash
rm *.so
/usr/bin/python2.7 setup.py build_ext


mv build/lib.linux-x86_64-2.7/df.so .
