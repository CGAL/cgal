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
           -Iinclude \
           -I../../include \
           -I../../../OpenNL/include \
           $(CGAL_CXXFLAGS) \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
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

# General targets:

all: \
          polyhedron_ex_parameterization$(EXE_EXT) \
          Simple_parameterization$(EXE_EXT) \
          Authalic_parameterization$(EXE_EXT) \
          Square_border_parameterization$(EXE_EXT) \
          Taucs_parameterization$(EXE_EXT) \
          Mesh_cutting_parameterization$(EXE_EXT) \
          Complete_parameterization_example$(EXE_EXT)

clean: \
          polyhedron_ex_parameterization.clean \
          Simple_parameterization.clean \
          Authalic_parameterization.clean \
          Square_border_parameterization.clean \
          Taucs_parameterization.clean \
          Mesh_cutting_parameterization.clean \
          Complete_parameterization_example.clean

# polyhedron_ex_parameterization:

polyhedron_ex_parameterization_objs = 	\
          Mesh_cutter$(OBJ_EXT) \
          polyhedron_ex_parameterization$(OBJ_EXT)

polyhedron_ex_parameterization$(EXE_EXT): $(polyhedron_ex_parameterization_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_ex_parameterization$(EXE_EXT) $(polyhedron_ex_parameterization_objs) $(LDFLAGS)

polyhedron_ex_parameterization.clean:
	rm -f polyhedron_ex_parameterization$(EXE_EXT) $(polyhedron_ex_parameterization_objs)

# Tutorial:

Simple_parameterization$(EXE_EXT): Simple_parameterization$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Simple_parameterization Simple_parameterization$(OBJ_EXT) $(LDFLAGS)

Authalic_parameterization$(EXE_EXT): Authalic_parameterization$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Authalic_parameterization Authalic_parameterization$(OBJ_EXT) $(LDFLAGS)

Square_border_parameterization$(EXE_EXT): Square_border_parameterization$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Square_border_parameterization Square_border_parameterization$(OBJ_EXT) $(LDFLAGS)

Taucs_parameterization$(EXE_EXT): Taucs_parameterization$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Taucs_parameterization Taucs_parameterization$(OBJ_EXT) $(LDFLAGS)

Mesh_cutting_parameterization$(EXE_EXT): Mesh_cutting_parameterization$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Mesh_cutting_parameterization Mesh_cutting_parameterization$(OBJ_EXT) $(LDFLAGS)

Complete_parameterization_example$(EXE_EXT): Complete_parameterization_example$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Complete_parameterization_example Complete_parameterization_example$(OBJ_EXT) $(LDFLAGS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

