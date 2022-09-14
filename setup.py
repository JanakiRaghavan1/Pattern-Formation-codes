from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
#from sage.env import SAGE_LOCAL, SAGE_ROOT
import os
SAGE_LOCAL = os.environ["SAGE_LOCAL"]
SAGE_ROOT  = os.environ["SAGE_ROOT"]

# get the annotated file as well
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

extra_compile_args = ['-Ofast', '-Wstrict-aliasing', '-fopenmp'] #['-w'
extra_link_args    = ['-L' + SAGE_LOCAL + '/lib', '-fopenmp']

ext_modules = [
        Extension('FHN',
            ['FHN.pyx'],
            libraries=['gsl', 'cblas', 'atlas'],
            include_dirs=[numpy.get_include(),
            # '%s/include/csage'%SAGE_LOCAL,
            # '%s/include'%SAGE_LOCAL,
            '%s/include/python2.7'%SAGE_LOCAL,
            # '%s/lib/python/site-packages/numpy/core/include'%SAGE_LOCAL,
            # '%s/devel/sage/sage/ext'%SAGE_ROOT,
            '%s/devel/sage'%SAGE_ROOT,
            '%s/devel/sage/sage/gsl'%SAGE_ROOT],
            library_dirs=[SAGE_LOCAL + '/lib/'],
            extra_compile_args = extra_compile_args,
            extra_link_args = extra_link_args),
        ]

setup(
        cmdclass = {'build_ext': build_ext},
        ext_modules = ext_modules,
        include_dirs = [#'%s/include/csage'%SAGE_LOCAL,
        # '%s/include'%SAGE_LOCAL,
        # '%s/include/python2.7'%SAGE_LOCAL,
        # '%s/lib/python/site-packages/numpy/core/include'%SAGE_LOCAL,
        # '%s/devel/sage/sage/ext'%SAGE_ROOT,
        # '%s/devel/sage'%SAGE_ROOT,
        '%s/devel/sage/sage/gsl'%SAGE_ROOT,
        '']
)
