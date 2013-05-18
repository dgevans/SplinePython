from distutils.core import setup

#setup.py
from distutils.extension import Extension

import numpy

files = ["SplineEigen.cpp"]
# define the name of the extension to use
extension_name='Spline_cpp'
extension_version ='1.0'
# define the directories to search for include files
# to get this to work, you may need to include the path
# to your boost installation. Mine was in
# '/usr/local/include', hence the corresponding entry.
include_dirs = ['/opt/local/include', '.',numpy.get_include()]
library_dirs = ['/usr/lib','/opt/local/lib']
# define the libraries to link with the boost python library
libraries = ['boost_python']
SuiteSparseLibraries = ['cxsparse','cholmod','amd','camd','ccolamd','SuiteSparse','suitesparseconfig']
# define the source files for the extension
source_files = ['SplineEigen.cpp']
# create the extension and add it to the python distribution
setup(name='SplinePython',
      version=extension_version,
      packages=['Spline'],
      package_dir={'Spline': 'SplineLib'},
      ext_package="Spline",
      ext_modules=[Extension(extension_name, source_files, include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries )])
#setup(
#    name='SplinePython',
#    version='',
#    packages=['Spline'],
#    url='',
#    license='',
#    author='dgevans',
#    author_email='',
#    description='',
#    package_dir={'Spline': 'SplineLib'},
#    package_data={'Spline': ['Spline_cpp.so']},
#)
