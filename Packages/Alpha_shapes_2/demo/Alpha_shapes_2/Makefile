# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#CGAL_MAKEFILE = $(UTIL)/CGAL/CGAL-2.3/make/makefile_sparc_SunOS-5.6_g++-2.95.3_LEDA
include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#----------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local

CXXFLAGS =   -I ../../include \
	     -I ../../../Triangulation_2/include \
             $(CGAL_CXXFLAGS) 

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



