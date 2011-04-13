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
           $(TESTSUITE_CXXFLAGS) \
           $(CGAL_CXXFLAGS) \
           -g

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                all_furthest_neighbors_2_example \
                extremal_polygon_2_example \
                rectangular_p_center_2_example \
                sorted_matrix_search_example 

all_furthest_neighbors_2_example: all_furthest_neighbors_2_example.o
	$(CGAL_CXX) -o all_furthest_neighbors_2_example all_furthest_neighbors_2_example.o $(LDFLAGS)

extremal_polygon_2_example: extremal_polygon_2_example.o
	$(CGAL_CXX) -o extremal_polygon_2_example extremal_polygon_2_example.o $(LDFLAGS)

rectangular_p_center_2_example: rectangular_p_center_2_example.o
	$(CGAL_CXX) -o rectangular_p_center_2_example rectangular_p_center_2_example.o $(LDFLAGS)

sorted_matrix_search_example: sorted_matrix_search_example.o
	$(CGAL_CXX) -o sorted_matrix_search_example sorted_matrix_search_example.o $(LDFLAGS)

clean:
	rm -f *.o \
                   all_furthest_neighbors_2_example \
                   extremal_polygon_2_example \
                   rectangular_p_center_2_example \
                   sorted_matrix_search_example \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

