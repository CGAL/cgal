# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = $(UTIL)/CGAL/CGAL-I/make/makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66_LEDA

include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local

CXXFLAGS = $(CGAL_CXXFLAGS) \
	   -mcpu=ultrasparc -mtune=ultrasparc -DNDEBUG \
           -B /u/alcor/0/prisme_util/Binutils/SunOS/bin/

#-mcpu=ultrasparc -mtune=ultrasparc -O2 -DNDEBUG
#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
LDFLAGS = \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#


clean:
	/bin/rm -f *.o \
                   demo_alpha demo_weight \
                   core

demo_alpha: demo_alpha.o
	$(CGAL_CXX) -o demo_alpha demo_alpha.o $(CGAL_WINDOW_LDFLAGS)

demo_weight: demo_weight.o
	$(CGAL_CXX) -o demo_weight demo_weight.o $(CGAL_WINDOW_LDFLAGS)


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c  $<



