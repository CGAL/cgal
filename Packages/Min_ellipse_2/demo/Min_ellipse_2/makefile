# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(CGAL_CXXFLAGS) \
           -g

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                demo_Min_circle_2 \
                demo_Min_ellipse_2 

demo_Min_circle_2: demo_Min_circle_2.o
	$(CGAL_CXX) -o demo_Min_circle_2 demo_Min_circle_2.o $(LDFLAGS)

demo_Min_ellipse_2: demo_Min_ellipse_2.o
	$(CGAL_CXX) -o demo_Min_ellipse_2 demo_Min_ellipse_2.o $(LDFLAGS)

clean:
	/bin/rm -f *.o \
                   demo_Min_circle_2 \
                   demo_Min_ellipse_2 \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

