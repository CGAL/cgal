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
                all_furthest_neighbors_2_test \
                extremal_polygon_2_test \
                rectangular_p_center_2_random1_test \
                rectangular_p_center_2_random2_test \
                sorted_matrix_search_test 

all_furthest_neighbors_2_test: all_furthest_neighbors_2_test.o
	$(CGAL_CXX) -o all_furthest_neighbors_2_test all_furthest_neighbors_2_test.o $(LDFLAGS)

extremal_polygon_2_test: extremal_polygon_2_test.o
	$(CGAL_CXX) -o extremal_polygon_2_test extremal_polygon_2_test.o $(LDFLAGS)

rectangular_p_center_2_random1_test: rectangular_p_center_2_random1_test.o
	$(CGAL_CXX) -o rectangular_p_center_2_random1_test rectangular_p_center_2_random1_test.o $(LDFLAGS)

rectangular_p_center_2_random2_test: rectangular_p_center_2_random2_test.o
	$(CGAL_CXX) -o rectangular_p_center_2_random2_test rectangular_p_center_2_random2_test.o $(LDFLAGS)

sorted_matrix_search_test: sorted_matrix_search_test.o
	$(CGAL_CXX) -o sorted_matrix_search_test sorted_matrix_search_test.o $(LDFLAGS)

clean:
	rm -f *.o \
                   all_furthest_neighbors_2_test \
                   extremal_polygon_2_test \
                   rectangular_p_center_2_random1_test \
                   rectangular_p_center_2_random2_test \
                   sorted_matrix_search_test \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

