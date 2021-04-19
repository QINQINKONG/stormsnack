from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
	"wetbulb",
        ["wetbulb.pyx"],
        libraries=["m"],
        extra_compile_args=['-fopenmp','-Ofast'],
        extra_link_args=['-fopenmp'],
    )
]

setup(
    name='wetbulb',
    ext_modules=cythonize(ext_modules,annotate=True),
    include_dirs=[numpy.get_include()],
    zip_safe=False,
)
