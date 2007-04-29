# Created by the script cgal_create_makefile and modified!!!!!
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#--------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)

PROFOPT= -g -fpermissive
CFLAGS=${PROFOPT} 
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
INTROSPEC_DIR = ../../../introspect
INTROSPEC_INC=-I$(INTROSPEC_DIR)

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   ${CFLAGS}  $(INTROSPEC_INC)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH) $(GLUT_LIBS_PATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

INTRO_LD_FLAGS= -L$(INTROSPEC_DIR) -lIntrospect -lglut

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:	visu

visu: 	SketchSample$(OBJ_EXT) visu$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ridges_viewer$(EXE_EXT) \
	SketchSample$(OBJ_EXT) visu$(OBJ_EXT) \
	$(LDFLAGS) $(INTRO_LD_FLAGS)

clean:  visu.clean SketchSample.clean 

depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# DO NOT DELETE

SketchSample.o: visu_poly.h enriched_polyhedron.h
visu_poly.o: enriched_polyhedron.h
