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

CXXFLAGS =  \
	-I../include/CGAL/IO \
        -I../include \
	-I../src/ \
	-I$(HOME)/Cgal/Viewer_3/include \
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
##---------------------------------------------------------------------#

all :    triangulation cycle_carre kernel_objects polyhedron



cycle_carre:  cycle_carre.o  
	$(CGAL_CXX) -DDEBUG -o cycle_carre cycle_carre.o  $(CGAL_WINDOW_LDFLAGS)

triangulation:  ./triangulation.o
	$(CGAL_CXX) -g -o ./triangulation ./triangulation.o $(VIEWER3_LDFLAGS)

kernel_objects:	./kernel_objects.o
	$(CGAL_CXX) -g -o ./kernel_objects ./kernel_objects.o $(VIEWER3_LDFLAGS)

polyhedron:  ./polyhedron.o
	$(CGAL_CXX) -g -o ./polyhedron ./polyhedron.o $(VIEWER3_LDFLAGS)


clean:
	/bin/rm -f *.o *.*~ *.ps triangulation cycle_carre kernel_objects polyhedron\
                core 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -g -c  $<




