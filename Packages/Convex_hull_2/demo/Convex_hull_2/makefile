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
                ch_of_polyline \
                convex_hull_2_demo \
                convex_hull_of_ipoints 

ch_of_polyline: ch_of_polyline.o
	$(CGAL_CXX) -o ch_of_polyline ch_of_polyline.o $(LDFLAGS)

convex_hull_2_demo: convex_hull_2_demo.o
	$(CGAL_CXX) -o convex_hull_2_demo convex_hull_2_demo.o $(LDFLAGS)

convex_hull_of_ipoints: convex_hull_of_ipoints.o
	$(CGAL_CXX) -o convex_hull_of_ipoints convex_hull_of_ipoints.o $(LDFLAGS)

clean: \
       ch_of_polyline.clean \
       convex_hull_2_demo.clean \
       convex_hull_of_ipoints.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

