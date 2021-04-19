from setuptools import setup
from Cython.Build import cythonize

setup(
    name='integrate',
    ext_modules=cythonize("integrate.pyx"),
    zip_safe=False,
)