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
           -DNDEBUG \
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
                test_Circular_kernel_arrangement$(EXE_EXT) \
                test_Circular_kernel_arrangement_Line_arc$(EXE_EXT) \
                test_Circular_kernel_arrangement_Line_Circular$(EXE_EXT) \
                test_Lazy_Circular_kernel_arrangement_Line_Circular$(EXE_EXT) \
                test_Lazy_Circular_kernel_arrangement_Line_arc$(EXE_EXT) \
                test_Lazy_Circular_kernel_arrangement$(EXE_EXT) \
                test_Circular_kernel$(EXE_EXT) \
                test_Line_arc$(EXE_EXT) 

test_Circular_kernel_arrangement$(EXE_EXT): test_Circular_kernel_arrangement$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Circular_kernel_arrangement test_Circular_kernel_arrangement$(OBJ_EXT) $(LDFLAGS)

test_Circular_kernel_arrangement_Line_arc$(EXE_EXT): test_Circular_kernel_arrangement_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Circular_kernel_arrangement_Line_arc test_Circular_kernel_arrangement_Line_arc$(OBJ_EXT) $(LDFLAGS)

test_Circular_kernel_arrangement_Line_Circular$(EXE_EXT): test_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Circular_kernel_arrangement_Line_Circular test_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT) $(LDFLAGS)

test_Lazy_Circular_kernel_arrangement_Line_Circular$(EXE_EXT): test_Lazy_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Lazy_Circular_kernel_arrangement_Line_Circular test_Lazy_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT) $(LDFLAGS)

test_Lazy_Circular_kernel_arrangement$(EXE_EXT): test_Lazy_Circular_kernel_arrangement$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Lazy_Circular_kernel_arrangement test_Lazy_Circular_kernel_arrangement$(OBJ_EXT) $(LDFLAGS)

test_Lazy_Circular_kernel_arrangement_Line_arc$(EXE_EXT): test_Lazy_Circular_kernel_arrangement_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Lazy_Circular_kernel_arrangement_Line_arc test_Lazy_Circular_kernel_arrangement_Line_arc$(OBJ_EXT) $(LDFLAGS)

test_Circular_kernel$(EXE_EXT): test_Circular_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Circular_kernel test_Circular_kernel$(OBJ_EXT) $(LDFLAGS)

test_Line_arc$(EXE_EXT): test_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Line_arc test_Line_arc$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_Circular_kernel_arrangement.clean \
                   test_Circular_kernel_arrangement_Line_arc.clean \
                   test_Circular_kernel_arrangement_Line_Circular.clean \
                   test_Lazy_Circular_kernel_arrangement_Line_Circular.clean \
                   test_Lazy_Circular_kernel_arrangement_Line_arc.clean \
                   test_Lazy_Circular_kernel_arrangement.clean \
                   test_Circular_kernel.clean \
                   test_Line_arc.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

