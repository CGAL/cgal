# This is the makefile for compiling a CGAL application.  
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#CGAL_MAKEFILE = /u/alcor/0/prisme_util/CGAL/CGAL-last/make/makefile_sparc_SunOS-5.6_g++-2.8.1_LEDA
#CGAL_MAKEFILE = $(UTIL)/CGAL/CGAL-last/make/makefile_sparc_SunOS-5.6_eg++-egcs-2.91.60_LEDA
#CGAL_MAKEFILE = /u/alcor/0/prisme_util/CGAL/CGAL-last/make/makefile_sparc_SunOS-5.6_CC-4.2_LEDA

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	-I../include \
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

test_triangulation_3: test_triangulation_3.o 
	$(CGAL_CXX) -o test_triangulation_3 test_triangulation_3.o $(LDFLAGS)

test_triangulation_geom_traits : test_triangulation_geom_traits.o
	$(CGAL_CXX)  -o test_triangulation_geom_traits test_triangulation_geom_traits.o $(LDFLAGS)

test_triangulation_tds: test_triangulation_tds.o
	$(CGAL_CXX)  -o test_triangulation_tds test_triangulation_tds.o $(LDFLAGS)

test_delaunay_3: test_delaunay_3.o 
	$(CGAL_CXX) -o test_delaunay_3 test_delaunay_3.o $(LDFLAGS)

clean:
	- rm -f *.o test_triangulation_3 test_triangulation_geom_traits test_triangulation_tds test_delaunay_3


realclean: clean
	- rm -f Program* Error* Compiler* error* Test_triangulation_IO_3

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX)  $(CXXFLAGS) -c $<
