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
           -I../../include \
           -I../../../Algebraic_kernel/include \
           -I../../../synaps \
           -I../../../Arrangement_2/include \
	   -I../../../Interval_arithmetic/include \
	   -I../../../Kernel_23/include \
	   -g \
	   $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
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
                test_Line_arc$(EXE_EXT) \
		test_Filtered_hexagon_kernel$(EXE_EXT) \
		test_Lazy_curved_kernel$(EXE_EXT) \
		test_Lazy_Line_arc$(EXE_EXT)

test_Circular_kernel$(EXE_EXT): test_Circular_kernel$(OBJ_EXT)
	$(CGAL_CXX) -pg $(LIBPATH) $(EXE_OPT)test_Circular_kernel test_Circular_kernel$(OBJ_EXT) $(LDFLAGS)

test_Filtered_hexagon_kernel$(EXE_EXT): test_Filtered_hexagon_kernel$(OBJ_EXT)
	$(CGAL_CXX) -pg $(LIBPATH) $(EXE_OPT)test_Filtered_hexagon_kernel test_Filtered_hexagon_kernel$(OBJ_EXT) $(LDFLAGS)

test_Lazy_curved_kernel$(EXE_EXT): test_Lazy_curved_kernel$(OBJ_EXT)
	$(CGAL_CXX) -pg $(LIBPATH) $(EXE_OPT)test_Lazy_curved_kernel test_Lazy_curved_kernel$(OBJ_EXT) $(LDFLAGS)

test_Line_arc$(EXE_EXT): test_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) -pg $(LIBPATH) $(EXE_OPT)test_Line_arc test_Line_arc$(OBJ_EXT) $(LDFLAGS)

test_Lazy_Line_arc$(EXE_EXT): test_Lazy_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) -pg $(LIBPATH) $(EXE_OPT)test_Lazy_Line_arc test_Lazy_Line_arc$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_Circular_kernel.clean \
		   test_Filtered_hexagon_kernel.clean \
                   test_Line_arc.clean \
		   test_Lazy_curved_kernel.clean \
		   test_Lazy_Line_arc.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

