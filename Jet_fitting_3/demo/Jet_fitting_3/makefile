# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#MARC
#CGAL_MAKEFILE = ${CGAL}/make/makefile_i686_Linux-2.6.12-10-386_g++-4.0.2
#CGAL_MAKEFILE = /home/mpouget/libs/CGAL-3.1/make/makefile_i686_Linux-2.6.12-10-386_g++-4.0.2

include $(CGAL_MAKEFILE)

PROFOPT= -g -fpermissive
CFLAGS=${PROFOPT} 
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
GSL_DIR=/usr/local/gsl
#GLUT_DIR=/user/sloriot/home/C/lib/glut-3.7
#MARC works without GLUT_DIR!! just needs GLUT_LIBS=-lGL -lGLU -lglut  -lX11 -lXmu -lXi -lm
#GLUT_DIR=/usr/local/glut-3.7
#GLUT_DIR=/home/mpouget/Desktop/glut-3.7

GSL_INC=-I$(GSL_DIR)/include

INTROSPEC_DIR = ../../../introspect
INTROSPEC_INC=-I$(INTROSPEC_DIR)

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   ${GSL_INC} $(GLUT_INC)\
           ${CFLAGS}  $(INTROSPEC_INC)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#


#LIBSX11= -L/usr/X11R6/lib
#include ${CGAL}/make/makefile_i686_Linux-2.4_g++-3.3.3 
##-H
##-D_VERBOSE_
#PROFOPT=-g -w
#CFLAGS=${PROFOPT} 
#ALL_INCLUDES=${GLUT_INC}


#GLUT_LIBS_PATH=-L$(GLUT_DIR)/lib/glut -L/usr/X11R6/lib
GLUT_LIBS=-lGL -lGLU -lglut  -lX11 -lXmu -lXi -lm
GLUT_LIBS=

GSL_LIBS_PATH=-L${GSL_DIR}/lib
GSL_LIBS=-lgsl -lgslcblas -lm

LIBPATH = \
          $(CGAL_LIBPATH) $(GSL_LIBS_PATH) $(GLUT_LIBS_PATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
          ${GSL_LIBS}  $(GLUT_LIBS)


EXE_EXT = .exe

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#


#VISU_OBJS=ogl_visu.o visu.o
#visu: $(VISU_OBJS)
#	$(CGAL_CXX) $(LIBPATH) -o visu.exe $(VISU_OBJS) \
#	$(LDFLAGS)

### visu with Camille/Laurent's tool
#test: test_visu.o
#	$(CGAL_CXX) $(LIBPATH) test_visu.o -o test $(LDFLAGS)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) -lGL -lglut

QTLDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) -lGL -lglut


INTRO_LD_FLAGS=-L$(INTROSPEC_DIR) -lInstrospect

all:            \
                visu


visu: SketchSample$(OBJ_EXT) visu$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)introspect-qt$(EXE_EXT) $(OBJECTS) SketchSample$(OBJ_EXT) visu$(OBJ_EXT) \
	$(QTLDFLAGS) $(INTRO_LD_FLAGS)

clean: 
	rm *$(EXE_EXT) *$(OBJ_EXT)

depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# DO NOT DELETE

SketchSample.o: SketchSample.h visu_poly.h enriched_polyhedron.h
SketchSample.o: visu_poly.h enriched_polyhedron.h
visu.o: SketchSample.h visu_poly.h enriched_polyhedron.h
visu_poly.o: enriched_polyhedron.h
