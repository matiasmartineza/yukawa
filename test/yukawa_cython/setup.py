from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize
import os
import numpy as np

ext = Extension("P2Pwrap", sources=["P2Pwrap.pyx"], language = "c", include_dirs = [np.get_include()], extra_compile_args=['-fopenmp'], extra_link_args=['-fopenmp'],)
setup(name="P2Pwrap", ext_modules = cythonize([ext]))

ext = Extension("calculateMultipoles", sources=["calculateMultipoles.pyx"], language = "c", include_dirs = [np.get_include()], extra_compile_args=['-fopenmp'], extra_link_args=['-fopenmp'],)
setup(name="calculateMultipoles", ext_modules = cythonize([ext]))

ext = Extension("M2Pwrap", sources=["M2Pwrap.pyx"], language = "c", include_dirs = [np.get_include()], extra_compile_args=['-fopenmp'], extra_link_args=['-fopenmp'],)
setup(name="M2Pwrap", ext_modules = cythonize([ext]))
