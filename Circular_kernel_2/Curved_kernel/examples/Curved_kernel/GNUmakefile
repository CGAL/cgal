# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
CGAL_MAKEFILE = /user/teillaud/home/CGAL/install/CGAL-I/make/makefile_i686_Linux-2.6_g++-3.4.4
#makefile_i686_Linux-2.6_g++33-3.3.4
#makefile_i686_Linux-2.6_g++-3.4.4
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include \
           -I../../../Algebraic_kernel/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
	example_arrangement_random_circles$(EXE_EXT) \
	example_arrangement_random_circles_segments$(EXE_EXT) \
	example_Exact_circular_kernel$(EXE_EXT) \
        example_Bbox_Hexagon_Lazy_Circular_kernel_arrangement_Line_Circular$(EXE_EXT) \
        example_Circular_kernel_arrangement$(EXE_EXT) \
        example_Circular_kernel_arrangement_Line_arc$(EXE_EXT) \
        example_Circular_kernel_arrangement_Line_Circular$(EXE_EXT) \
        example_Lazy_Circular_kernel_arrangement$(EXE_EXT) \
        example_Lazy_Circular_kernel_arrangement_Line_arc$(EXE_EXT) \
        example_Lazy_Circular_kernel_arrangement_Line_Circular$(EXE_EXT) 

example_arrangement_random_circles$(EXE_EXT): example_arrangement_random_circles$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_arrangement_random_circles example_arrangement_random_circles$(OBJ_EXT) $(LDFLAGS)

example_arrangement_random_circles_segments$(EXE_EXT): example_arrangement_random_circles_segments$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_arrangement_random_circles_segments example_arrangement_random_circles_segments$(OBJ_EXT) $(LDFLAGS)

example_Exact_circular_kernel$(EXE_EXT): example_Exact_circular_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Exact_circular_kernel example_Exact_circular_kernel$(OBJ_EXT) $(LDFLAGS)

example_Bbox_Hexagon_Lazy_Circular_kernel_arrangement_Line_Circular$(EXE_EXT): example_Bbox_Hexagon_Lazy_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Bbox_Hexagon_Lazy_Circular_kernel_arrangement_Line_Circular example_Bbox_Hexagon_Lazy_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT) $(LDFLAGS)

example_Circular_kernel_arrangement$(EXE_EXT): example_Circular_kernel_arrangement$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Circular_kernel_arrangement example_Circular_kernel_arrangement$(OBJ_EXT) $(LDFLAGS)

example_Circular_kernel_arrangement_Line_arc$(EXE_EXT): example_Circular_kernel_arrangement_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Circular_kernel_arrangement_Line_arc example_Circular_kernel_arrangement_Line_arc$(OBJ_EXT) $(LDFLAGS)

example_Circular_kernel_arrangement_Line_Circular$(EXE_EXT): example_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Circular_kernel_arrangement_Line_Circular example_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT) $(LDFLAGS)

example_Lazy_Circular_kernel_arrangement$(EXE_EXT): example_Lazy_Circular_kernel_arrangement$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Lazy_Circular_kernel_arrangement example_Lazy_Circular_kernel_arrangement$(OBJ_EXT) $(LDFLAGS)

example_Lazy_Circular_kernel_arrangement_Line_arc$(EXE_EXT): example_Lazy_Circular_kernel_arrangement_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Lazy_Circular_kernel_arrangement_Line_arc example_Lazy_Circular_kernel_arrangement_Line_arc$(OBJ_EXT) $(LDFLAGS)

example_Lazy_Circular_kernel_arrangement_Line_Circular$(EXE_EXT): example_Lazy_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_Lazy_Circular_kernel_arrangement_Line_Circular example_Lazy_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT) $(LDFLAGS)

clean: \
	example_arrangement_random_circles.clean \
        example_arrangement_random_circles_segments.clean \
	example_Exact_circular_kernel.clean \
	example_Bbox_Hexagon_Lazy_Circular_kernel_arrangement_Line_Circular.clean \
        example_Circular_kernel_arrangement.clean \
        example_Circular_kernel_arrangement_Line_arc.clean \
        example_Circular_kernel_arrangement_Line_Circular.clean \
        example_Lazy_Circular_kernel_arrangement.clean \
        example_Lazy_Circular_kernel_arrangement_Line_arc.clean \
        example_Lazy_Circular_kernel_arrangement_Line_Circular.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

