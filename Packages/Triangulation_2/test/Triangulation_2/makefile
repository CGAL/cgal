# This is the makefile for compiling a CGAL application.  
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = /u/alcor/0/prisme_util/CGAL/CGAL-1.2-I-4/make/makefile_sparc_SunOS-5.6_g++-2.8.1_LEDA
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	-I../../include \
	-I./include \
	$(CGAL_CXXFLAGS)\
	$(TESTSUITE_CXXFLAGS) \
	-g


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: test_triangulation_geom_traits test_triangulation_tds \
	test_triangulation_2 test_delaunay_triangulation_2 test_constrained_triangulation_2

test_triangulation_geom_traits: test_triangulation_geom_traits.o
	$(CGAL_CXX) -o test_triangulation_geom_traits test_triangulation_geom_traits.o $(LDFLAGS)

test_triangulation_tds: test_triangulation_tds.o
	$(CGAL_CXX) -o test_triangulation_tds test_triangulation_tds.o $(LDFLAGS)

test_triangulation_2: test_triangulation_2.o
	$(CGAL_CXX) -o test_triangulation_2 test_triangulation_2.o $(LDFLAGS)

test_delaunay_triangulation_2: test_delaunay_triangulation_2.o
	$(CGAL_CXX) -o test_delaunay_triangulation_2 test_delaunay_triangulation_2.o $(LDFLAGS)

test_constrained_triangulation_2: test_constrained_triangulation_2.o	
	$(CGAL_CXX) -o test_constrained_triangulation_2 test_constrained_triangulation_2.o $(LDFLAGS)

clean:
	- rm -f *.o test_triangulation_geom_traits test_triangulation_2 \
	            test_triangulation_tds test_delaunay_triangulation_2 core \
	            T??.triangulation test_constrained_triangulation_2

realclean: clean
	- rm -f Program* Error* Compiler* error*

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX)  $(CXXFLAGS) -c $<
