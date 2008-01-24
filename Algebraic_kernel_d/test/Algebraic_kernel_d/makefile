# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
CGAL_CXXFLAGS = -ftemplate-depth-100 -Wall -frounding-math -ffriend-injection -ffriend-injection -I/KM/home-0/emeliyan/work/CGAL-3.4-I-69/include/CGAL/config/i686_Linux-2.6_g++-4.1.2 -I/KM/home-0/emeliyan/work/CGAL-3.4-I-69/include -I/KM/home-0/emeliyan/gcc/include/c++/4.1.2 -pthread -I/KM/home-0/emeliyan/gcc/include/c++/4.1.2 -pthread -It -I/usr/X11R6/include -I/KM/home-0/emeliyan/thirdparty_libs/gmp/include -I/KM/home-0/emeliyan/thirdparty_libs/gmp/include -I/KM/home-0/emeliyan/thirdparty_libs/mpfr/include -I/LEDA/LEDA-5.x/LEDA-5.2-std/linux/g++-4.1/incl -I/KM/projects/ecg/software/linux/g++-4.1/qt-x11-free-3.3.6/include

CGAL_LDFLAGS = -L/KM/home-0/emeliyan/work/CGAL-3.4-I-69/lib/i686_Linux-2.6_g++-4.1.2 -L/usr/X11R6/lib -L/KM/home-0/emeliyan/thirdparty_libs/gmp/lib -L/KM/home-0/emeliyan/thirdparty_libs/gmp/lib -L/KM/home-0/emeliyan/thirdparty_libs/mpfr/lib -L/LEDA/LEDA-5.x/LEDA-5.2-std/linux/g++-4.1 -L/KM/projects/ecg/software/linux/g++-4.1/qt-x11-free-3.3.6/lib -Wl,-R/KM/home-0/emeliyan/work/CGAL-3.4-I-69/lib/i686_Linux-2.6_g++-4.1.2:/usr/X11R6/lib:/KM/home-0/emeliyan/thirdparty_libs/gmp/lib:/KM/home-0/emeliyan/thirdparty_libs/gmp/lib:/KM/home-0/emeliyan/thirdparty_libs/mpfr/lib:/LEDA/LEDA-5.x/LEDA-5.2-std/linux/g++-4.1:/KM/projects/ecg/software/linux/g++-4.1/qt-x11-free-3.3.6/lib -L/KM/home-0/emeliyan/gcc/lib -lpthread -lCGAL -lCGALcore++ -lCGALQt -lCGALimageIO -lCGALPDB -lqt-mt -lGL -lGLU -lP -lG -lL -lz -lmpfr -lgmpxx -lgmp -lX11 -lm

CXXFLAGS = \
           -Iinclude \
           -I../../include \
	   -I/KM/home-0/emeliyan/tmp/CGAL_packages/Installation/include \
	   -I/KM/home-0/emeliyan/tmp/CGAL_packages/STL_Extension/include \
	   -I/KM/home-0/emeliyan/tmp/CGAL_packages/Modular_arithmetic/include \
	   -I/KM/home-0/emeliyan/tmp/CGAL_packages/Polynomial/include \
	   -I/KM/home-0/emeliyan/tmp/CGAL_packages/Algebraic_foundations/include \
	   -I/KM/home-0/emeliyan/tmp/CGAL_packages/Number_types/include \
	   -I/KM/home-0/emeliyan/tmp/CGAL_packages/Number_types/test/Number_types/include \
	   -I/KM/home-0/emeliyan/tmp/EXACUS/NumeriX/include \
	   -I/KM/home-0/emeliyan/tmp/EXACUS/SweepX/include \
	   -I/KM/home-0/emeliyan/tmp/EXACUS/AlciX/include \
	   -I/KM/home-0/emeliyan/tmp/EXACUS/Support/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = -L/tmp/emeliyan/exacus/lib -lNiX -Wl,--rpath -Wl,/tmp/emeliyan/exacus/lib  \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            simple_models$(EXE_EXT)  

Algebraic_curve_kernel_2$(EXE_EXT): Algebraic_curve_kernel_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Algebraic_curve_kernel_2 Algebraic_curve_kernel_2$(OBJ_EXT) $(LDFLAGS)

clean: Algebraic_curve_kernel_2.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

