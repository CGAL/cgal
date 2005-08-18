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

XXFLAGS = \
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

all:      \
          polyhedron_ex_parameterization$(EXE_EXT) \
          Polyhedron_parameterization1$(EXE_EXT) \
          Polyhedron_parameterization2$(EXE_EXT) \
          Polyhedron_parameterization3$(EXE_EXT) \
          Polyhedron_parameterization4$(EXE_EXT) \
          Polyhedron_parameterization5$(EXE_EXT) \
          Polyhedron_parameterization6$(EXE_EXT)

clean:    \
          polyhedron_ex_parameterization.clean \
          Polyhedron_parameterization1.clean \
          Polyhedron_parameterization2.clean \
          Polyhedron_parameterization3.clean \
          Polyhedron_parameterization4.clean \
          Polyhedron_parameterization5.clean \
          Polyhedron_parameterization6.clean


# polyhedron_ex_parameterization:

polyhedron_ex_parameterization_objs = 	\
          Mesh_cutter$(OBJ_EXT) \
          options$(OBJ_EXT) \
          polyhedron_ex_parameterization$(OBJ_EXT)

polyhedron_ex_parameterization$(EXE_EXT): $(polyhedron_ex_parameterization_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_ex_parameterization$(EXE_EXT) $(polyhedron_ex_parameterization_objs) $(LDFLAGS)

polyhedron_ex_parameterization.clean:
	rm polyhedron_ex_parameterization$(EXE_EXT) $(polyhedron_ex_parameterization_objs)


# Polyhedron_parameterization1:

Polyhedron_parameterization1_objs = 	\
          Polyhedron_parameterization1$(OBJ_EXT)

Polyhedron_parameterization1$(EXE_EXT): $(Polyhedron_parameterization1_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polyhedron_parameterization1$(EXE_EXT) $(Polyhedron_parameterization1_objs) $(LDFLAGS)

Polyhedron_parameterization1.clean:
	rm Polyhedron_parameterization1$(EXE_EXT) $(Polyhedron_parameterization1_objs)


# Polyhedron_parameterization2:

Polyhedron_parameterization2_objs = 	\
          Polyhedron_parameterization2$(OBJ_EXT)

Polyhedron_parameterization2$(EXE_EXT): $(Polyhedron_parameterization2_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polyhedron_parameterization2$(EXE_EXT) $(Polyhedron_parameterization2_objs) $(LDFLAGS)

Polyhedron_parameterization2.clean:
	rm Polyhedron_parameterization2$(EXE_EXT) $(Polyhedron_parameterization2_objs)


# Polyhedron_parameterization3:

Polyhedron_parameterization3_objs = 	\
          Polyhedron_parameterization3$(OBJ_EXT)

Polyhedron_parameterization3$(EXE_EXT): $(Polyhedron_parameterization3_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polyhedron_parameterization3$(EXE_EXT) $(Polyhedron_parameterization3_objs) $(LDFLAGS)

Polyhedron_parameterization3.clean:
	rm Polyhedron_parameterization3$(EXE_EXT) $(Polyhedron_parameterization3_objs)


# Polyhedron_parameterization4:

Polyhedron_parameterization4_objs = 	\
          Polyhedron_parameterization4$(OBJ_EXT)

Polyhedron_parameterization4$(EXE_EXT): $(Polyhedron_parameterization4_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polyhedron_parameterization4$(EXE_EXT) $(Polyhedron_parameterization4_objs) $(LDFLAGS)

Polyhedron_parameterization4.clean:
	rm Polyhedron_parameterization4$(EXE_EXT) $(Polyhedron_parameterization4_objs)


# Polyhedron_parameterization5:

Polyhedron_parameterization5_objs = 	\
          Polyhedron_parameterization5$(OBJ_EXT)

Polyhedron_parameterization5$(EXE_EXT): $(Polyhedron_parameterization5_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polyhedron_parameterization5$(EXE_EXT) $(Polyhedron_parameterization5_objs) $(LDFLAGS)

Polyhedron_parameterization5.clean:
	rm Polyhedron_parameterization5$(EXE_EXT) $(Polyhedron_parameterization5_objs)


# Polyhedron_parameterization6:

Polyhedron_parameterization6_objs = 	\
          Polyhedron_parameterization6$(OBJ_EXT)

Polyhedron_parameterization6$(EXE_EXT): $(Polyhedron_parameterization6_objs)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polyhedron_parameterization6$(EXE_EXT) $(Polyhedron_parameterization6_objs) $(LDFLAGS)

Polyhedron_parameterization6.clean:
	rm Polyhedron_parameterization6$(EXE_EXT) $(Polyhedron_parameterization6_objs)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

