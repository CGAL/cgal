

FLTK_DIR = ${UTIL}/FLTK/Linux/fltk-1.0.6
GL_DIR = /0/prisme_util/MESA/Mesa-3.1-Linux

CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-I/make/makefile_i686_Linux-2.2.14_g++-2.95.2_LEDA
include $(CGAL_MAKEFILE)



#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CHECKS_OFF = -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS \
	    -DCGAL_NO_ASSERTIONS -DCGAL_NO_WARNINGS


CXXFLAGS = -g  \
	   -I${FLTK_DIR} \
	   -I${GL_DIR}/include \
           -I../../include \
	   -I/usr/X11R6/include \
	   -B$(GCC_EXEC_PREFIX) $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(CGAL_CXXFLAGS) $(CHECKS_OFF) -DUSE_THREAD

#-I${GL_DIR}/include \

LDFLAGS =  $(CGAL_LDFLAGS) -L${GL_DIR}/lib -L$(FLTK_DIR)/lib -L/usr/X11R6/lib \
		-Wl,-R$(FLTK_DIR)/lib -Wl,-R${GL_DIR}/lib -Wl,-R/usr/X11R6/lib \
	 	-lfltk -lXext -lX11 -lGL -lGLU

#-lMesaGL -lMesaGLU

#-L$(GL_DIR)/lib 
#-R$(GL_DIR)/lib 
#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#


demo: main.o 
	$(CGAL_CXX) -o $@ $(OFILES) main.o -lpthread $(LDFLAGS)

facet: facet_obj.o
	$(CGAL_CXX) -o $@ $(OFILES) facet_obj.o -lpthread $(LDFLAGS)


tetra1: tetra1.o
	$(CGAL_CXX) -o $@ $(OFILES) tetra1.o -lpthread $(LDFLAGS)

tetra2: tetra2.o
	$(CGAL_CXX) -o $@ $(OFILES) tetra2.o -lpthread $(LDFLAGS)

tetra3: tetra3.o
	$(CGAL_CXX) -o $@ $(OFILES) tetra3.o -lpthread $(LDFLAGS)

terrain: terrain.o
	$(CGAL_CXX) -o $@ $(OFILES) terrain.o -lpthread $(LDFLAGS)

all: facet demo tetra1 tetra2 tetra3

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c  -D_REENTRANT $<

.c.o:
	$(CGAL_CXX) $(CXXFLAGS) -c -D_REENTRANT  $<




clean:
	rm -f *.o core *~ 
	rm -f ii_files/*
distclean: clean
	rm -f demo facet tetra1 tetra2 tetra3 
	rm ../include/CGAL/*~
