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

# Turn on optimizer and extensive checking
EXTRA_FLAGS = \
           -O3 \
           -DCGAL_PARAMETERIZATION_CHECK_EXPENSIVE

CXXFLAGS = \
           -I../../include \
           -I../../../OpenNL/include \
           $(CGAL_CXXFLAGS) \
           $(TESTSUITE_CXXFLAGS) \
           -I$(TAUCSROOT)/src \
           -I$(TAUCSROOT)/build/$(TAUCS_OSTYPE) \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH) \
          -L$(TAUCSROOT)/lib/$(TAUCS_OSTYPE) \
          -L$(TAUCSROOT)/external/lib/$(TAUCS_OSTYPE)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
          -ltaucs -llapack -lf77blas -lcblas -latlas -lmetis -lg2c


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:      polyhedron_ex_parameterization$(EXE_EXT)

clean:	  polyhedron_ex_parameterization.clean

polyhedron_ex_parameterization_objs = 	\
          Mesh_cutter$(OBJ_EXT) \
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

