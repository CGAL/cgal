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
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           -g

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
           $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                in_place_list_prog \
                min_element_if_example_noheader \
                min_max_element_example_noheader 

in_place_list_prog: in_place_list_prog.o
	$(CGAL_CXX) -o in_place_list_prog in_place_list_prog.o $(LDFLAGS)

min_element_if_example_noheader: min_element_if_example_noheader.o
	$(CGAL_CXX) -o min_element_if_example_noheader min_element_if_example_noheader.o $(LDFLAGS)

min_max_element_example_noheader: min_max_element_example_noheader.o
	$(CGAL_CXX) -o min_max_element_example_noheader min_max_element_example_noheader.o $(LDFLAGS)

clean:
	rm -f *.o \
                   in_place_list_prog \
                   min_element_if_example_noheader \
                   min_max_element_example_noheader \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

