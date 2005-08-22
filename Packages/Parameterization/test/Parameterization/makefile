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
           -Iinclude \
           -I../../examples/Parameterization \
           -I../../examples/Parameterization/include \
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

all:      \
          extensive_parameterization_test$(EXE_EXT)

clean:    \
          extensive_parameterization_test.clean

# extensive_parameterization_test:

extensive_parameterization_test_objs = 	\
          Mesh_cutter$(OBJ_EXT) \
          extensive_parameterization_test$(OBJ_EXT)

#Mesh_cutter$(OBJ_EXT) : ../../examples/Parameterization/Mesh_cutter.C
Mesh_cutter$(OBJ_EXT) :
	echo '*****'
	pwd
	echo '*****'
	ls -l ..
	echo '*****'
	ls -l ../..
	echo '*****'
	echo Content of CGAL_INCL_DIR/../examples is:
	ls $(CGAL_INCL_DIR)/../examples
	echo '*****'
	env | sort
	echo '*****'
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) ../../examples/Parameterization/Mesh_cutter.C

extensive_parameterization_test$(EXE_EXT): $(extensive_parameterization_test_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)extensive_parameterization_test$(EXE_EXT) $(extensive_parameterization_test_objs) $(LDFLAGS)

extensive_parameterization_test.clean:
	rm -f extensive_parameterization_test$(EXE_EXT) $(extensive_parameterization_test_objs)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

