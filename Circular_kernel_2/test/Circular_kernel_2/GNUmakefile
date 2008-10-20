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
           -I../../../Intersections_2/include \
           -I../../include \
           -Iinclude\
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
                test_Circular_kernel$(EXE_EXT) \
		test_Exact_circular_kernel$(EXE_EXT)\
		test_Filtered_bbox_circular_kernel$(EXE_EXT)\
                test_Line_arc$(EXE_EXT)

test_Circular_kernel$(EXE_EXT): test_Circular_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Circular_kernel test_Circular_kernel$(OBJ_EXT) $(LDFLAGS)

test_Exact_circular_kernel$(EXE_EXT): test_Exact_circular_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Exact_circular_kernel test_Exact_circular_kernel$(OBJ_EXT) $(LDFLAGS)

test_Filtered_bbox_circular_kernel$(EXE_EXT): test_Filtered_bbox_circular_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Filtered_bbox_circular_kernel test_Filtered_bbox_circular_kernel$(OBJ_EXT) $(LDFLAGS)

test_Line_arc$(EXE_EXT): test_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Line_arc test_Line_arc$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_Circular_kernel.clean \
		   test_Exact_circular_kernel.clean \
		   test_Filtered_bbox_circular_kernel.clean\
                   test_Line_arc.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
