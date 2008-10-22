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

CXXFLAGS = \
           -I../../../Kernel_23/include \
           -I../../../Cartesian_kernel/include \
           -I../../../Homogeneous_kernel/include \
           -I../../../Intersections_3/include \
           -I../../include \
           -Iinclude\
           -I../../../Kernel_23/test/Kernel_23/include \
           -I../Kernel_23/include \
           -I../Cartesian_kernel/include \
           -I../Homogeneous_kernel/include \
           -I../Intersections_3/include \
           -g\
	         $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

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
                test_Spherical_kernel$(EXE_EXT) \
		test_Exact_spherical_kernel$(EXE_EXT)\
		test_Lazy_Spherical_kernel$(EXE_EXT)\
		test_Spherical_kernel_basics$(EXE_EXT)
				
test_Spherical_kernel$(EXE_EXT): test_Spherical_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Spherical_kernel test_Spherical_kernel$(OBJ_EXT) $(LDFLAGS)

test_Exact_spherical_kernel$(EXE_EXT): test_Exact_spherical_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Exact_spherical_kernel test_Exact_spherical_kernel$(OBJ_EXT) $(LDFLAGS)

test_Lazy_Spherical_kernel$(EXE_EXT): test_Lazy_Spherical_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Lazy_Spherical_kernel test_Lazy_Spherical_kernel$(OBJ_EXT) $(LDFLAGS)

test_Spherical_kernel_basics$(EXE_EXT): test_Spherical_kernel_basics$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Spherical_kernel_basics test_Spherical_kernel_basics$(OBJ_EXT) $(LDFLAGS)

clean: \
       test_Spherical_kernel.clean \
       test_Spherical_kernel_basics.clean \
		   test_Exact_spherical_kernel.clean \
		   test_Lazy_Spherical_kernel.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
