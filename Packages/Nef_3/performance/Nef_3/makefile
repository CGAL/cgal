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

#           -O2 -DNDEBUG \

CXXFLAGS = \
	   $(ADDITIONAL_CXXFLAGS) \
           -DCGAL_NEF3_DUMP_STATISTICS \
           -DCGAL_NEF3_TIMER_INITIALIZE_KDTREE \
           -DCGAL_NEF3_TIMER_PLANE_SWEEPS \
           -DCGAL_NEF3_TIMER_SPHERE_SWEEPS \
           -DCGAL_NEF3_TIMER_POINT_LOCATION \
           -DCGAL_NEF3_TIMER_SIMPLIFICATION \
           -DCGAL_NEF3_TIMER_INTERSECTION \
           -DCGAL_NEF3_TIMER_EXTERNAL_STRUCTURE \
           -DCGAL_NEF3_TIMER_OVERLAY \
	   -DCGAL_NEF3_TIMER_BINARY_OPERATION \
           -I/home/hachenb/CGAL/Nef_3/include \
           -I/home/hachenb/CGAL/Nef_S2/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(DEBUG_OPT)

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
                union_of_rotated_cylinders$(EXE_EXT) \
	        two_tetrahedra_fields$(EXE_EXT) \
	       	tetrahedra_and_grid$(EXE_EXT) \
                big_and_small$(EXE_EXT) \
                quadratic$(EXE_EXT)

union_of_rotated_cylinders$(EXE_EXT): union_of_rotated_cylinders$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)union_of_rotated_cylinders union_of_rotated_cylinders$(OBJ_EXT) $(LDFLAGS)
two_tetrahedra_fields$(EXE_EXT): two_tetrahedra_fields$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)two_tetrahedra_fields two_tetrahedra_fields$(OBJ_EXT) $(LDFLAGS)
tetrahedra_and_grid$(EXE_EXT): tetrahedra_and_grid$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tetrahedra_and_grid tetrahedra_and_grid$(OBJ_EXT) $(LDFLAGS)
big_and_small$(EXE_EXT): big_and_small$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)big_and_small big_and_small$(OBJ_EXT) $(LDFLAGS)
quadratic$(EXE_EXT): quadratic$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)quadratic quadratic$(OBJ_EXT) $(LDFLAGS)

clean: \
                   union_of_rotated_cylinders.clean \
	           two_tetrahedra_fields.clean \
	           tetrahedra_and_grid.clean \
	           big_and_small.clean \
	           quadratic.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

