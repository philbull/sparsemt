from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
    name = "sparse",
    version = "0.1",
    description = "ghost analytic galaxy-halo model",
    ext_modules = [ 
                Extension(
                          "sparse",
                          ["sparse.c", ],
                          libraries = ['m', 'gsl', 'gslcblas', 'gomp'],
                          depends = ['sparse.h'],
                          extra_compile_args=['-O4', '-fopenmp', '-std=c99']
                          ) ],
    include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs(),
)
