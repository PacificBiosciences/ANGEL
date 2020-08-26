from setuptools import Extension, setup
from Cython.Build import cythonize

ext_modules = [Extension("c_ORFscores", sources=["c_ORFscores.pyx"], language="c++")]

setup(
		name = 'c_ORF',
		ext_modules = cythonize(ext_modules)
)

