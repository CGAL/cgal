
CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-136/make/makefile_i686_CYGWINNT-5.0-1.3.12_CL.EXE-1300

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
# -DCGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS

CXXFLAGS = \
	   -Id:/cgal/develop/Intersections_2/include \
	   -Id:/cgal/develop/Cartesian_kernel/include \
	   -I../../include \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
test_delaunay_2 \
test_polyhedron_prog_incr_builder \
MyKernel



MyKernel$(EXE_EXT): MyKernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)MyKernel MyKernel$(OBJ_EXT) $(LDFLAGS)

test_delaunay_2$(EXE_EXT): test_delaunay_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_delaunay_2 test_delaunay_2$(OBJ_EXT) $(LDFLAGS)

test_polyhedron_prog_incr_builder$(EXE_EXT): test_polyhedron_prog_incr_builder$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_polyhedron_prog_incr_builder test_polyhedron_prog_incr_builder$(OBJ_EXT) $(LDFLAGS)

clean: \
MyKernel.clean \
test_delaunay_2.clean \
test_polyhedron_prog_incr_builder.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

