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
                all_furthest_neighbors_2_demo \
                extremal_polygon_2_demo \
                rectangular_p_center_2_demo 

all_furthest_neighbors_2_demo: all_furthest_neighbors_2_demo.o
	$(CGAL_CXX) -o all_furthest_neighbors_2_demo all_furthest_neighbors_2_demo.o $(LDFLAGS)

extremal_polygon_2_demo: extremal_polygon_2_demo.o
	$(CGAL_CXX) -o extremal_polygon_2_demo extremal_polygon_2_demo.o $(LDFLAGS)

rectangular_p_center_2_demo: rectangular_p_center_2_demo.o
	$(CGAL_CXX) -o rectangular_p_center_2_demo rectangular_p_center_2_demo.o $(LDFLAGS)

clean:
	rm -f *.o \
                   all_furthest_neighbors_2_demo \
                   extremal_polygon_2_demo \
                   rectangular_p_center_2_demo \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

