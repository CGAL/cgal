# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#CGAL_MAKEFILE=$(UTIL)/CGAL/CGAL-I/make/makefile_sparc_SunOS-5.6_g++-2.95.2_LEDA
CGAL_MAKEFILE = $(UTIL)/CGAL/CGAL-2.2/make/makefile_sparc_SunOS-5.6_g++-2.95.2_LEDA
#makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66_LEDA
#makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66_LEDA
#makefile_i686_Linux-2.2.14_g++-2.95.2_LEDA
#makefile_sparc_SunOS-5.6_g++-2.95.2_LEDA
#makefile_sparc_SunOS-5.6_g++-2.95.2_LEDA

include $(CGAL_MAKEFILE)

FLTK_DIR = ${UTIL}/FLTK/SunOS/fltk-1.0.6
# ---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local

CHECKS_OFF = -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS \
	    -DCGAL_NO_ASSERTIONS -DCGAL_NO_WARNINGS

CXXFLAGS =  \
	-I../../include \
	-I../../../Postscript/include/CGAL/IO \
	-I../../../Postscript/src/ \
	-I../include/CGAL/IO \
	-I../include \
	-I../src/ \
	-I${FLTK_DIR} \
	-I/usr/X11R6/include/ \
	-B$(UTIL)/Binutils/$(UNAME)/bin/ \
	-B$(GCC_EXEC_PREFIX) $(LONG_NAME_PROBLEM_CXXFLAGS) \
	$(CGAL_CXXFLAGS) $(CHECKS_OFF) -DUSE_THREAD


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

VIEWER3_LDFLAGS = \
	$(CGAL_WINDOW_LDFLAGS) \
	 -L$(FLTK_DIR)/lib -L/usr/lib/X11 \
	-R$(FLTK_DIR)/lib -R/usr/lib/X11 -Wl \
	 -lfltk -lXext -lX11 -lGL -lGLU -lpthread


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#


demo: demo.o 
	$(CGAL_CXX) -o $@ $(OFILES) demo.o -lpthread $(LDFLAGS)


facet_obj: facet_obj.o
	$(CGAL_CXX) -o $@ $(OFILES) facet_obj.o -lpthread $(LDFLAGS)


tetra1: tetra1.o
	$(CGAL_CXX) -o $@ $(OFILES) tetra1.o -lpthread $(LDFLAGS)

tetra2: tetra2.o
	$(CGAL_CXX) -o $@ $(OFILES) tetra2.o -lpthread $(LDFLAGS)

tetra3: tetra3.o
	$(CGAL_CXX) -o $@ $(OFILES) tetra3.o -lpthread $(LDFLAGS)

terrain: terrain.o
	$(CGAL_CXX) -o $@ $(OFILES) terrain.o -lpthread $(LDFLAGS)

testter: testter.o
	$(CGAL_CXX) -o $@ $(OFILES) testter.o -lpthread $(LDFLAGS)

stream_ex: stream_ex.o
	$(CGAL_CXX) -o $@ $(OFILES) stream_ex.o -lpthread $(LDFLAGS)

all: facet_obj demo tetra1 tetra2 tetra3 terrain stream_ex

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
	rm -f demo facet tetra1 tetra2 tetra3 terrain stream_ex
	rm ../include/CGAL/*~
