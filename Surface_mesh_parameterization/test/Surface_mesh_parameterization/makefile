# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    location of other Surface_mesh_parameterization subfolders
#                    (non standard when running the test suite)
#---------------------------------------------------------------------#

CGAL = $(CGAL_INCL_DIR)/..

EXAMPLES = `./find_first ../../examples/Surface_mesh_parameterization $(CGAL)/examples/Surface_mesh_parameterization`

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -Iinclude \
           -I$(EXAMPLES) \
           -I$(EXAMPLES)/include \
           -I../../include \
           -I../../../OpenNL/include \
           -I../../../Installation/include \
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

all:      \
          extensive_parameterization_test$(EXE_EXT)

clean:    \
          extensive_parameterization_test.clean

# extensive_parameterization_test:

Mesh_cutter$(OBJ_EXT) :
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(EXAMPLES)/Mesh_cutter.cpp

extensive_parameterization_test_objs = 	\
          Mesh_cutter$(OBJ_EXT) \
          extensive_parameterization_test$(OBJ_EXT)

extensive_parameterization_test$(EXE_EXT): $(extensive_parameterization_test_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)extensive_parameterization_test$(EXE_EXT) $(extensive_parameterization_test_objs) $(LDFLAGS)

extensive_parameterization_test.clean:
	rm -f extensive_parameterization_test$(EXE_EXT) $(extensive_parameterization_test_objs)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

