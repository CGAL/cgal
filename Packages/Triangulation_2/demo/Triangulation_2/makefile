# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)



#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I./include 	-I../../include \
	$(CGAL_CXXFLAGS) \
	$(LONG_NAME_PROBLEM_CXXFLAGS) \
	-DDELAUNAY


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(CGAL_LDFLAGS) \
	$(LONG_NAME_PROBLEM_LDFLAGS) 
WINDOW_LDFLAGS = $(CGAL_WINDOW_LDFLAGS) 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#



demo: demo.o parse.o
	$(CGAL_CXX) -o demo demo.o parse.o $(WINDOW_LDFLAGS)



clean:
	- rm -rf *~ *.o core tr

realclean: clean
	- rm -rf lfc demo grid crust points3D voronoi \
	triangulation_basic_example triangulation_terrain \
	triangulation_locate_grid

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
