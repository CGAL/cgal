# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
INTROSPEC_DIR = ../../../introspect
INTROSPEC_INC=-I$(INTROSPEC_DIR)

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(INTROSPEC_INC)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
#GLUT_LIBS_PATH=-L$(GLUT_DIR)/lib/glut -L/usr/X11R6/lib
GLUT_LIBS=-lGL -lGLU -lglut  -lX11 -lXmu -lXi -lm

LIBPATH = \
          $(CGAL_LIBPATH) $(GLUT_LIBS_PATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
          $(GLUT_LIBS)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) -lGL -lglut

QTLDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) -lGL -lglut


INTRO_LD_FLAGS=-L$(INTROSPEC_DIR) -lIntrospect

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:	visu

visu: 	SketchSample$(OBJ_EXT) visu$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)introspect-qt \
	SketchSample$(OBJ_EXT) visu$(OBJ_EXT) \
	$(QTLDFLAGS) $(INTRO_LD_FLAGS)

clean: 	SketchSample.clean \
	visu.clean

depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# DO NOT DELETE

SketchSample.o: SketchSample.h visu_poly.h enriched_polyhedron.h
SketchSample.o: visu_poly.h enriched_polyhedron.h
visu.o: SketchSample.h visu_poly.h enriched_polyhedron.h
visu_poly.o: enriched_polyhedron.h
