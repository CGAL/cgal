# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
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
#GLUT_LIBS_PATH=-L$(GLUT_DIR)/lib/glut -L/usr/X11R6/lib
GLUT_LIBS=-lGL -lGLU -lglut  -lX11 -lXmu -lXi -lm

LIBPATH = \
          $(CGAL_LIBPATH) $(GLUT_LIBS_PATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
          $(GLUT_LIBS)

EXE_EXT = .exe

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) -lGL -lglut

QTLDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) -lGL -lglut


INTRO_LD_FLAGS=-L$(INTROSPEC_DIR) -lInstrospect


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:	visu

visu: 	SketchSample$(OBJ_EXT) visu$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)introspect-qt$(EXE_EXT) \
	SketchSample$(OBJ_EXT) visu$(OBJ_EXT) \
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

ogl_visu.o: ogl_visu.h /usr/include/math.h /usr/include/features.h
ogl_visu.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
ogl_visu.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
ogl_visu.o: /usr/include/bits/mathcalls.h /usr/include/GL/glut.h
ogl_visu.o: /usr/include/GL/gl.h /usr/include/GL/glext.h
ogl_visu.o: /usr/include/GL/glu.h enriched_polyhedron.h
ogl_visu.o: /usr/include/math.h /usr/include/features.h
ogl_visu.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
ogl_visu.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
ogl_visu.o: /usr/include/bits/mathcalls.h /usr/include/GL/glut.h
ogl_visu.o: /usr/include/GL/gl.h /usr/include/GL/glext.h
ogl_visu.o: /usr/include/GL/glu.h enriched_polyhedron.h
SketchSample.o: SketchSample.h visu_poly.h enriched_polyhedron.h
SketchSample.o: visu_poly.h enriched_polyhedron.h
visu.o: SketchSample.h visu_poly.h enriched_polyhedron.h
visu_poly.o: enriched_polyhedron.h
