# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right inclcoude file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE

CGAL_MAKEFILE = /usr/local/CGAL/CGAL-2.3-I/make/makefile_i686_Linux-2.2.14-5.0smp_g++-2.95-2.95.2_LEDA

include $(CGAL_MAKEFILE)
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	$(CGAL_CXXFLAGS) \
	$(LONG_NAME_PROBLEM_CXXFLAGS) \
	$(DEBUG_OPT) \
	-g

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
	$(CGAL_LIBPATH)

LDFLAGS = \
	$(LONG_NAME_PROBLEM_LDFLAGS) \
	$(CGAL_LDFLAGS) \
	$(CGAL_WINDOW_LDFLAGS)\
	$(TESTSUITE_LDFLAGS) \

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:	EmptyRect

EmptyRect: EmptyRect.o
	$(CGAL_CXX) $(CXXFLAGS) $(LIBPATH) -o EmptyRect EmptyRect.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f EmptyRect
	rm -f core
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<










