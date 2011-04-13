# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)
#include /KM/projects/CGAL/Work_in_Progress/CGAL-2.1-I-32/make/makefile_sparc_SunOS-5.6_g++-2.95.2_LEDA++
#include /KM/projects/CGAL/Work_in_Progress/CGAL-2.1-I-39/make/makefile_sparc_SunOS-5.6_g++-2.95.2_LEDA++

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(DEBUG_OPT)

#           -DCGAL_NO_LEDA_HANDLE \

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                test_kernel__

test_kernel__$(EXE_EXT): test_kernel__$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernel__ test_kernel__$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_kernel__.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

