# python setup.py build_ext --inplace

from distutils.core import setup
from Cython.Build import cythonize

setup(name='myPairwiseAlignment',
      ext_modules=cythonize("myPairwiseAlignment.pyx"))