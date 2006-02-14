

FLTK_DIR = ${UTIL}/FLTK/SunOS/fltk-1.0.6
#GL_DIR = /0/prisme_util/MESA/SunOS/Mesa-3.0/include/
# GL_DIR = /usr/local/mesa
CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-I/make/makefile_sparc_SunOS-5.6_g++-2.95.2_LEDA
include $(CGAL_MAKEFILE)



#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CHECKS_OFF = -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS \
	    -DCGAL_NO_ASSERTIONS -DCGAL_NO_WARNINGS


CXXFLAGS = -g  \
	   -I${FLTK_DIR} \
	   -I/usr/X11R6/include/ \
	   -B$(UTIL)/Binutils/SunOS/bin/ \
	   -B$(GCC_EXEC_PREFIX) $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(CGAL_CXXFLAGS) $(CHECKS_OFF) -DCGAL_USE_FLTK -DCGAL_USE_OPENGL

#-I${GL_DIR}/include \

LDFLAGS =  $(CGAL_LDFLAGS) -L$(FLTK_DIR)/lib -L/usr/lib/X11 \
		-R$(FLTK_DIR)/lib -R/usr/lib/X11 \
	 	-lfltk -lXext -lX11 -lGL -lGLU

#-lMesaGL -lMesaGLU

#-L$(GL_DIR)/lib 
#-R$(GL_DIR)/lib 
#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#


test: Viewer_test.o 
	$(CGAL_CXX) -o $@ $(OFILES) Viewer_test.o  $(LDFLAGS)

all: test

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c   $<

.c.o:
	$(CGAL_CXX) $(CXXFLAGS) -c   $<




clean:
	rm -f *.o core *~
	rm -f ii_files/*
	rm -f test error.txt
	rm -f CompilerOutput*
	rm -f ErrorOutput_sparc_SunOS-5.6_g++-2.95.2_LEDA



