# Created by the script create_makefile
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

CXXFLAGS = \
           -I../../../_2/include \
           -I../../../_3/include \
           -I../../../H3/include \
           -I../../../H2/include \
           -I../../../Cartesian_kernel/include \
           -I../../../Kernel_basic/include \
           -I../../../Number_types/include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                test_kernel__$(EXE_EXT) 

test_kernel__$(EXE_EXT): test_kernel__$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernel__ test_kernel__$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_kernel__.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

