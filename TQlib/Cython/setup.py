from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
#This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

#compile the fortran modules without linking
fortran_mod_comp = 'gfortran -I../../ ../C/liboctqisoc.F90 -c -o liboctqisoc.o -O3 -fPIC'
print fortran_mod_comp
system(fortran_mod_comp)
shared_obj_comp = 'gfortran -I../../ ../liboctq.F90 -c -o liboctq.o -O3 -fPIC'
print shared_obj_comp
system(shared_obj_comp)

ext_modules=[Extension(# module name:
			'pyoctq',
			#source file:
			['pyoctq.pyx'],
			#other compile args for gcc
			extra_compile_args=['-fPIC', '-O3'],
			#other files to link to
			extra_link_args=['liboceq.a', 'liboctqisoc.o', 'liboctq.o','-lgfortran'])]
setup(name = 'pyoctq',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with Numpy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)

