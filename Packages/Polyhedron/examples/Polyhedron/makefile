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
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
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
                polyhedron_prog_incr_builder \
                polyhedron_prog_normals \
                polyhedron_prog_off \
                polyhedron_prog_point_iterator \
                polyhedron_prog_simple \
                polyhedron_prog_tetra 

polyhedron_prog_incr_builder$(EXE_EXT): polyhedron_prog_incr_builder$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_incr_builder polyhedron_prog_incr_builder$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_normals$(EXE_EXT): polyhedron_prog_normals$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_normals polyhedron_prog_normals$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_off$(EXE_EXT): polyhedron_prog_off$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_off polyhedron_prog_off$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_point_iterator$(EXE_EXT): polyhedron_prog_point_iterator$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_point_iterator polyhedron_prog_point_iterator$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_simple$(EXE_EXT): polyhedron_prog_simple$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_simple polyhedron_prog_simple$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_tetra$(EXE_EXT): polyhedron_prog_tetra$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_tetra polyhedron_prog_tetra$(OBJ_EXT) $(LDFLAGS)

clean: \
                   polyhedron_prog_incr_builder.clean \
                   polyhedron_prog_normals.clean \
                   polyhedron_prog_off.clean \
                   polyhedron_prog_point_iterator.clean \
                   polyhedron_prog_simple.clean \
                   polyhedron_prog_tetra.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

