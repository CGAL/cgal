# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#FC
#CGAL_MAKEFILE = ${CGAL}/make/makefile_i686_Linux-2.6.11-1.35-FC3_g++-3.4.4
#MARC
#CGAL_MAKEFILE = ${CGAL}/make/makefile_i686_Linux-2.6.12-10-386_g++-4.0.2
CGAL_MAKEFILE = /home/mpouget/libs/CGAL-3.1/make/makefile_i686_Linux-2.6.12-10-386_g++-4.0.2

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
#GLUT_INC=-I$(GLUT_DIR)/include

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   ${GSL_INC} $(GLUT_INC)\
           ${CFLAGS} 

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

GSL_LIBS_PATH=-L${GSL_DIR}/lib
GSL_LIBS=-lgsl -lgslcblas -lm

LIBPATH = \
          $(CGAL_LIBPATH) $(GSL_LIBS_PATH) $(GLUT_LIBS_PATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
          ${GSL_LIBS}  $(GLUT_LIBS)


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                visu

VISU_OBJS=ogl_visu.o visu.o
visu: $(VISU_OBJS)
	$(CGAL_CXX) $(LIBPATH) -o visu.exe $(VISU_OBJS) \
	$(LDFLAGS)

rmexe:
	rm *.exe

clean: \
        visu.clean \
        ogl_visu.clean \
	rmexe


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
visu.o: ogl_visu.h /usr/include/math.h /usr/include/features.h
visu.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
visu.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
visu.o: /usr/include/bits/mathcalls.h /usr/include/GL/glut.h
visu.o: /usr/include/GL/gl.h /usr/include/GL/glext.h /usr/include/GL/glu.h
visu.o: enriched_polyhedron.h
