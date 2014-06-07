from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("c_ORFscores", ["c_ORFscores.pyx"], language="c++")]

setup(
		name = 'c_ORF',
		cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules
)

