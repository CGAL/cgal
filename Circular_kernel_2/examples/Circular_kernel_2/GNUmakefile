# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
#makefile_i686_Linux-2.6_g++33-3.3.4
#makefile_i686_Linux-2.6_g++-3.4.4
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include \
           -I../../../Algebraic_kernel_for_circles/include \
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
	arrangement_random_circles$(EXE_EXT) \
	arrangement_random_circles_segments$(EXE_EXT) \
	Exact_circular_kernel$(EXE_EXT) \
        Circular_kernel_arrangement$(EXE_EXT) \
        Circular_kernel_arrangement_Line_arc$(EXE_EXT) \
        Circular_kernel_arrangement_Line_Circular$(EXE_EXT) \
        intersecting_arcs$(EXE_EXT) \
        Lazy_Circular_kernel_arrangement$(EXE_EXT) \
        Lazy_Circular_kernel_arrangement_Line_arc$(EXE_EXT) \
        Lazy_Circular_kernel_arrangement_Line_Circular$(EXE_EXT) 

arrangement_random_circles$(EXE_EXT): arrangement_random_circles$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)arrangement_random_circles arrangement_random_circles$(OBJ_EXT) $(LDFLAGS)

arrangement_random_circles_segments$(EXE_EXT): arrangement_random_circles_segments$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)arrangement_random_circles_segments arrangement_random_circles_segments$(OBJ_EXT) $(LDFLAGS)

Exact_circular_kernel$(EXE_EXT): Exact_circular_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Exact_circular_kernel Exact_circular_kernel$(OBJ_EXT) $(LDFLAGS)

Circular_kernel_arrangement$(EXE_EXT): Circular_kernel_arrangement$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Circular_kernel_arrangement Circular_kernel_arrangement$(OBJ_EXT) $(LDFLAGS)

Circular_kernel_arrangement_Line_arc$(EXE_EXT): Circular_kernel_arrangement_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Circular_kernel_arrangement_Line_arc Circular_kernel_arrangement_Line_arc$(OBJ_EXT) $(LDFLAGS)

Circular_kernel_arrangement_Line_Circular$(EXE_EXT): Circular_kernel_arrangement_Line_Circular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Circular_kernel_arrangement_Line_Circular Circular_kernel_arrangement_Line_Circular$(OBJ_EXT) $(LDFLAGS)

intersecting_arcs$(EXE_EXT): intersecting_arcs$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)intersecting_arcs intersecting_arcs$(OBJ_EXT) $(LDFLAGS)

Lazy_Circular_kernel_arrangement$(EXE_EXT): Lazy_Circular_kernel_arrangement$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Lazy_Circular_kernel_arrangement Lazy_Circular_kernel_arrangement$(OBJ_EXT) $(LDFLAGS)

Lazy_Circular_kernel_arrangement_Line_arc$(EXE_EXT): Lazy_Circular_kernel_arrangement_Line_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Lazy_Circular_kernel_arrangement_Line_arc Lazy_Circular_kernel_arrangement_Line_arc$(OBJ_EXT) $(LDFLAGS)

Lazy_Circular_kernel_arrangement_Line_Circular$(EXE_EXT): Lazy_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Lazy_Circular_kernel_arrangement_Line_Circular Lazy_Circular_kernel_arrangement_Line_Circular$(OBJ_EXT) $(LDFLAGS)

clean: \
	arrangement_random_circles.clean \
        arrangement_random_circles_segments.clean \
	Exact_circular_kernel.clean \
        Circular_kernel_arrangement.clean \
        Circular_kernel_arrangement_Line_arc.clean \
        Circular_kernel_arrangement_Line_Circular.clean \
        intersecting_arcs.clean \
        Lazy_Circular_kernel_arrangement.clean \
        Lazy_Circular_kernel_arrangement_Line_arc.clean \
        Lazy_Circular_kernel_arrangement_Line_Circular.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

