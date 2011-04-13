# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = $(UTIL)/CGAL/CGAL-2.3-I-84/make/makefile_sparc_SunOS-5.6_g++-2.95.3_LEDA

#$(UTIL)/CGAL-I/make/makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66_LEDA
#makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66_LEDA
#makefile_sparc_SunOS-5.6_g++-2.95.1_LEDA
#makefile_sparc_SunOS-5.6_g++-2.95.2_LEDA
#makefile_i686_Linux-2.2.10_eg++-egcs-2.91.66_LEDA
#makefile_i686_Linux-2.2.10_g++-2.95.1_LEDA
#makefile_i686_Linux-2.2.10_g++-2.95.2_LEDA

include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#----------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local

CXXFLAGS =   -I$(HOME)/Cgal-I/For_Submission/Alpha_shapes_2/include \
             $(CGAL_CXXFLAGS) \
	     -mcpu=ultrasparc -mtune=ultrasparc -O2 \
	     -B /u/alcor/0/prisme_util/Binutils/SunOS/bin/

#-mcpu=ultrasparc -mtune=ultrasparc -Wall -W -O2 -DDEBUG
#-B /u/alcor/0/prisme_util/Binutils/SunOS/bin/ -g
#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
LDFLAGS = $(CGAL_WINDOW_LIBPATH) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: demo_alpha demo_weight

demo_alpha: demo_alpha.o
	$(CGAL_CXX) -o demo_alpha demo_alpha.o $(LDFLAGS)

demo_weight: demo_weight.o
	$(CGAL_CXX) -o demo_weight demo_weight.o $(LDFLAGS)

clean:
	\rm -f *.o  *.obj *.exe \
                   demo_alpha demo_weight \
                   core


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c  $<



