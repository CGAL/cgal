# This is the makefile for compiling a CGAL application.  
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.



#CGAL_MAKEFILE = /u/alcor/0/prisme_util/CGAL/CGAL-I/make/makefile_i686_Linux-2.2.7_eg++-egcs-2.91.66_LEDA
CGAL_MAKEFILE = /u/alcor/0/prisme_util/CGAL/CGAL-2.1-I-11/make/makefile_i686_Linux-2.2.7_g++-2.95._LEDA/
include $(CGAL_MAKEFILE)


#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	-I../../include \
	-I./include \
	$(TESTSUITE_CXXFLAGS) \
	$(CGAL_CXXFLAGS)\
	-DCGAL_NO_PRECONDITIONS -DCGAL_NO_ASSERTIONS -DCGAL_NO_POSTCONDITIONS
#	-O
#	-g 
#	-B$(UTIL)/Binutils/SunOS/bin/ 
#	$(LONG_NAME_PROBLEM_CXXFLAGS) 
#	-g \
#       -B$(UTIL)/Binutils/SunOS/bin/ 
#	-DCGAL_NO_PRECONDITIONS -DCGAL_NO_ASSERTIONS -DCGAL_NO_POSTCONDITIONS
#	-fhonor-std
#	-I /u/rigel/0/prisme/spion/prog/Cgal/Triangulation/include \


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)
#	  $(LONG_NAME_PROBLEM_LDFLAGS) \

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: test_triangulation_geom_traits test_triangulation_tds \
	test_triangulation_2 test_delaunay_triangulation_2 test_constrained_triangulation_2

test_triangulation_geom_traits: test_triangulation_geom_traits.o
	$(CGAL_CXX) -o test_triangulation_geom_traits test_triangulation_geom_traits.o $(LDFLAGS)

test_regular_triangulation_2: test_regular_triangulation_2.o
	$(CGAL_CXX) -o test_regular_triangulation_2 test_regular_triangulation_2.o $(LDFLAGS)

test_triangulation_tds: test_triangulation_tds.o
	$(CGAL_CXX) -o test_triangulation_tds test_triangulation_tds.o $(LDFLAGS)

test_triangulation_2: test_triangulation_2.o
	$(CGAL_CXX) -o test_triangulation_2 test_triangulation_2.o $(LDFLAGS)

test_terrain : test_terrain.o
	$(CGAL_CXX) -o test_terrain  test_terrain.o $(LDFLAGS)

terr_trian : terr_trian.o
	$(CGAL_CXX) -o terr_trian  terr_trian.o $(LDFLAGS)

essai : essai.o
	$(CGAL_CXX) -o essai  essai.o $(LDFLAGS)

test_delaunay_triangulation_2: test_delaunay_triangulation_2.o
	$(CGAL_CXX) -o test_delaunay_triangulation_2 test_delaunay_triangulation_2.o $(LDFLAGS)

test_constrained_triangulation_2: test_constrained_triangulation_2.o	
	$(CGAL_CXX) -o test_constrained_triangulation_2 test_constrained_triangulation_2.o $(LDFLAGS)

clean:
	- rm -f *.o test_triangulation_geom_traits test_triangulation_2 \
	            test_triangulation_tds test_delaunay_triangulation_2 core \
	            T??.triangulation  file_tds? \
		    test_constrained_triangulation_2 test_regular_triangulation_2

realclean: clean
	- rm -f Program* Error* Compiler* error*

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX)  $(CXXFLAGS) -c $<
