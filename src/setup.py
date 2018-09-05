from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext


#extensions = Extension("input1", sources=["input1.pyx", "initialize.c"])
ext_modules = [
    Extension("python_bridge", ["python_bridge.pyx"], include_dirs=['module/'])
    ]

setup(
  name = 'app',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules  
)