# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = $(UTIL)/CGAL/CGAL-2.0/make/makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66_LEDA

include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local

CXXFLAGS = -I../../include \
        -I./include \
	-I$(HOME)/Cgal/Postscript/include \
	$(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
LDFLAGS = \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

postscript: postscript.o  PS_Stream.o
	$(CGAL_CXX) -o postscript  PS_Stream.o postscript.o $(CGAL_WINDOW_LDFLAGS)

essai1: essai1.o
	$(CGAL_CXX) -o essai1  essai1.o $(CGAL_WINDOW_LDFLAGS)
	

clean:
	/bin/rm -f *.o \
                   postscript \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c  $<




