# This is the makefile for compiling a CGAL application.  
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.


CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-2.3-I-79/make/makefile_sparc_SunOS-5.6_g++-3.0.

include $(CGAL_MAKEFILE)


#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	-I./include \
	-I../../include \
	$(TESTSUITE_CXXFLAGS) \
	$(CGAL_CXXFLAGS)\
	-g 


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#


test_triangulation_geom_traits: test_triangulation_geom_traits.o
	$(CGAL_CXX) -o test_triangulation_geom_traits test_triangulation_geom_traits.o $(LDFLAGS)

test_triangulation_tds: test_triangulation_tds.o
	$(CGAL_CXX) -o test_triangulation_tds test_triangulation_tds.o $(LDFLAGS)

test_triangulation_tdsul: test_triangulation_tdsul.o
	$(CGAL_CXX) -o test_triangulation_tdsul test_triangulation_tdsul.o $(LDFLAGS)

test_triangulation_2: test_triangulation_2.o
	$(CGAL_CXX) -o test_triangulation_2 test_triangulation_2.o $(LDFLAGS)

test_triangulation_2_bis: test_triangulation_2_bis.o
	$(CGAL_CXX) -o test_triangulation_2_bis test_triangulation_2_bis.o $(LDFLAGS)

test_delaunay_triangulation_2: test_delaunay_triangulation_2.o
	$(CGAL_CXX) -o test_delaunay_triangulation_2 test_delaunay_triangulation_2.o $(LDFLAGS)

test_regular_triangulation_2: test_regular_triangulation_2.o
	$(CGAL_CXX) -o test_regular_triangulation_2 test_regular_triangulation_2.o $(LDFLAGS)

test_constrained_triangulation_2: test_constrained_triangulation_2.o	
	$(CGAL_CXX) -o test_constrained_triangulation_2 test_constrained_triangulation_2.o $(LDFLAGS)

test_const_del_triangulation_2: test_const_del_triangulation_2.o
	$(CGAL_CXX) -o test_const_del_triangulation_2 test_const_del_triangulation_2.o $(LDFLAGS)

clean:
	- rm -f *.o  *.obj *.exe \
	test_triangulation_geom_traits \
	test_triangulation_2 test_triangulation_2_bis\
	test_triangulation_tds test_triangulation_tdsul \
	test_triangulation_2_bis \
	test_delaunay_triangulation_2 \
	test_constrained_triangulation_2 \
	test_regular_triangulation_2 \
	test_const_del_triangulation_2 \
	 core \
	T??.triangulation  file_tds? 

realclean: clean
	- rm -f Program* Error* Compiler* error*

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX)  $(CXXFLAGS) -c $<
