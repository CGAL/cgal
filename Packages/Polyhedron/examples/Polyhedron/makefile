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
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                polyhedron_prog_color \
                polyhedron_prog_cube \
                polyhedron_prog_incr_builder \
                polyhedron_prog_incr_builder_old \
                polyhedron_prog_normals \
                polyhedron_prog_normals_old \
                polyhedron_prog_off \
                polyhedron_prog_off_old \
                polyhedron_prog_planes \
                polyhedron_prog_point_iterator_old \
                polyhedron_prog_simple \
                polyhedron_prog_simple_old \
                polyhedron_prog_subdiv \
                polyhedron_prog_subdiv_with_boundary \
                polyhedron_prog_tetra \
                polyhedron_prog_tetra_old \
                polyhedron_prog_traits \
                polyhedron_prog_vector 

polyhedron_prog_color$(EXE_EXT): polyhedron_prog_color$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_color polyhedron_prog_color$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_cube$(EXE_EXT): polyhedron_prog_cube$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_cube polyhedron_prog_cube$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_incr_builder$(EXE_EXT): polyhedron_prog_incr_builder$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_incr_builder polyhedron_prog_incr_builder$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_incr_builder_old$(EXE_EXT): polyhedron_prog_incr_builder_old$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_incr_builder_old polyhedron_prog_incr_builder_old$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_normals$(EXE_EXT): polyhedron_prog_normals$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_normals polyhedron_prog_normals$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_normals_old$(EXE_EXT): polyhedron_prog_normals_old$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_normals_old polyhedron_prog_normals_old$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_off$(EXE_EXT): polyhedron_prog_off$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_off polyhedron_prog_off$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_off_old$(EXE_EXT): polyhedron_prog_off_old$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_off_old polyhedron_prog_off_old$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_planes$(EXE_EXT): polyhedron_prog_planes$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_planes polyhedron_prog_planes$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_point_iterator_old$(EXE_EXT): polyhedron_prog_point_iterator_old$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_point_iterator_old polyhedron_prog_point_iterator_old$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_simple$(EXE_EXT): polyhedron_prog_simple$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_simple polyhedron_prog_simple$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_simple_old$(EXE_EXT): polyhedron_prog_simple_old$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_simple_old polyhedron_prog_simple_old$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_subdiv$(EXE_EXT): polyhedron_prog_subdiv$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_subdiv polyhedron_prog_subdiv$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_subdiv_with_boundary$(EXE_EXT): polyhedron_prog_subdiv_with_boundary$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_subdiv_with_boundary polyhedron_prog_subdiv_with_boundary$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_tetra$(EXE_EXT): polyhedron_prog_tetra$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_tetra polyhedron_prog_tetra$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_tetra_old$(EXE_EXT): polyhedron_prog_tetra_old$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_tetra_old polyhedron_prog_tetra_old$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_traits$(EXE_EXT): polyhedron_prog_traits$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_traits polyhedron_prog_traits$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_vector$(EXE_EXT): polyhedron_prog_vector$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_vector polyhedron_prog_vector$(OBJ_EXT) $(LDFLAGS)

clean: \
                   polyhedron_prog_color.clean \
                   polyhedron_prog_cube.clean \
                   polyhedron_prog_incr_builder.clean \
                   polyhedron_prog_incr_builder_old.clean \
                   polyhedron_prog_normals.clean \
                   polyhedron_prog_normals_old.clean \
                   polyhedron_prog_off.clean \
                   polyhedron_prog_off_old.clean \
                   polyhedron_prog_planes.clean \
                   polyhedron_prog_point_iterator_old.clean \
                   polyhedron_prog_simple.clean \
                   polyhedron_prog_simple_old.clean \
                   polyhedron_prog_subdiv.clean \
                   polyhedron_prog_subdiv_with_boundary.clean \
                   polyhedron_prog_tetra.clean \
                   polyhedron_prog_tetra_old.clean \
                   polyhedron_prog_traits.clean \
                   polyhedron_prog_vector.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

