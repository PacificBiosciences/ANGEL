import os
from setuptools import setup, find_packages
from distutils.extension import Extension

ext_modules = [Extension("Angel.c_ORFscores", ["src/c_ORFscores.cpp"], language="c++")]

setup(
        name = 'Angel',
        version = '1.0', 
        author = 'Elizabeth Tseng',
        author_email = 'etseng@pacificbiosciences.com',
        license = 'LICENSE.txt',
        ext_modules = ext_modules,
        scripts = ['dumb_predict.py', 'angel_train.py', 'angel_predict.py'], #TODO
        package_dir = {'Angel': 'src'}, 
        packages = ['Angel'],
        zip_safe = False,
        install_requires = [
            'numpy >= 1.7',
            'biopython >= 1.6',
            'scikit-learn',
            ]
)
