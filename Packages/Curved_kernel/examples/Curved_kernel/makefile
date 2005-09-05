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
                example_Circular_kernel_arrangement$(EXE_EXT) \
                example_Circular_kernel_arrangement_Line_arc$(EXE_EXT) \
                example_Circular_kernel_arrangement_Line_Circular$(EXE_EXT) \
                example_Lazy_Circular_kernel_arrangement_Line_Circular$(EXE_EXT) \
                example_Lazy_Circular_kernel_arrangement_Line_arc$(EXE_EXT) \
                example_Lazy_Circular_kernel_arrangement$(EXE_EXT)


example_Circular_kernel_arrangement$(EXE_EXT): example_Circular_kernel_arrangement$(OBJ_EXT)
	$(CGAL_CXX) -pg $(LIBPATH) $(EXE_OPT)example_Circular_kernel_arrangement example_Circular_kernel_arrangement$(OBJ_EXT) $(LDFLAGS)

example_Circular_kernel_arrangement_Line_arc$(EXE_EXT): example_Circular_kernel_arrangement_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) -pg $(LIBPATH) $(EXE_OPT)example_Circular_kernel_arrangement_Line_arc example_Circular_kernel_arrangement_Line_arc$(OBJ_EXT) $(LDFLAGS)

example_Circular_kernel_arrangement_Line_Circular$(EXE_EXT): example_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT)
	$(CGAL_CXX) -pg $(LIBPATH) $(EXE_OPT)example_Circular_kernel_arrangement_Line_Circular example_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT) $(LDFLAGS)

example_Lazy_Circular_kernel_arrangement_Line_Circular$(EXE_EXT): example_Lazy_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Lazy_Circular_kernel_arrangement_Line_Circular example_Lazy_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT) $(LDFLAGS)

example_Lazy_Circular_kernel_arrangement$(EXE_EXT): example_Lazy_Circular_kernel_arrangement$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Lazy_Circular_kernel_arrangement example_Lazy_Circular_kernel_arrangement$(OBJ_EXT) $(LDFLAGS)

example_Lazy_Circular_kernel_arrangement_Line_arc$(EXE_EXT): example_Lazy_Circular_kernel_arrangement_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Lazy_Circular_kernel_arrangement_Line_arc example_Lazy_Circular_kernel_arrangement_Line_arc$(OBJ_EXT) $(LDFLAGS)

clean: \
                   example_Circular_kernel_arrangement.clean \
                   example_Circular_kernel_arrangement_Line_arc.clean \
                   example_Circular_kernel_arrangement_Line_Circular.clean \
                   example_Lazy_Circular_kernel_arrangement_Line_Circular.clean \
                   example_Lazy_Circular_kernel_arrangement_Line_arc.clean \
                   example_Lazy_Circular_kernel_arrangement.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

