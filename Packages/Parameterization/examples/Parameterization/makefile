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
           -I../../../OpenNL/include \
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

all:            polyhedron_ex_parameterization$(EXE_EXT)

clean:			polyhedron_ex_parameterization.clean

polyhedron_ex_parameterization_objs = 	\
               	Mesh_cutter$(OBJ_EXT) \
                Mesh_feature_extractor$(OBJ_EXT) \
                options$(OBJ_EXT) \
                polyhedron_ex_parameterization$(OBJ_EXT) 

polyhedron_ex_parameterization$(EXE_EXT): $(polyhedron_ex_parameterization_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_ex_parameterization$(EXE_EXT) $(polyhedron_ex_parameterization_objs) $(LDFLAGS)

polyhedron_ex_parameterization.clean:
	rm polyhedron_ex_parameterization$(EXE_EXT) $(polyhedron_ex_parameterization_objs)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

