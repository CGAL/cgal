# This is the makefile for compiling a CGAL application.  
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#CGAL_MAKEFILE = /u/alcor/0/prisme_util/CGAL/CGAL-last/make/makefile_sparc_SunOS-5.6_g++-2.8.1_LEDA
#CGAL_MAKEFILE = $(UTIL)/CGAL/CGAL-last/make/makefile_sparc_SunOS-5.6_eg++-egcs-2.91.60_LEDA
#CGAL_MAKEFILE = /u/alcor/0/prisme_util/CGAL/CGAL-last/make/makefile_sparc_SunOS-5.6_CC-4.2_LEDA

# pour purify, virer LEDA
#CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-last/make/makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66

CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-2.1-I-17/make/makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66_LEDA

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

# pour purify :
#CGAL_CXX=purify /usr/local/egcs/bin/g++

CXXFLAGS = \
	-I../../include \
	-I./include \
        -B$(UTIL)/Binutils/$(UNAME)/bin/ \
	 $(CGAL_CXXFLAGS)\
	 $(TESTSUITE_CXXFLAGS) \
	-g
#         -gdwarf-2
#	-O2 -fno-default-inline \
#	-DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS -DCGAL_NO_ASSERTIONS
#         -g


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

# pour purify :
#CGAL_LDFLAGS= $(CGAL_LIB_DIR)/$(CGAL_OS_COMPILER)/libCGAL.a \
#              -L$(GMP_LIB_DIR) -lgmp -lm 

LDFLAGS = $(TESTSUITE_LDFLAGS) \
          -B$(UTIL)/Binutils/$(UNAME)/bin/ \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: test_triangulation_geom_traits test_triangulation_tds \
     test_triangulation_3 test_delaunay_3 test_regular_3 
	./test_triangulation_geom_traits
	./test_triangulation_tds
	./test_triangulation_3
	./test_delaunay_3
	./test_regular_3 

test_triangulation_geom_traits : test_triangulation_geom_traits.o
	$(CGAL_CXX)  -o test_triangulation_geom_traits test_triangulation_geom_traits.o $(LDFLAGS)

test_triangulation_tds: test_triangulation_tds.o
	$(CGAL_CXX)  -o test_triangulation_tds test_triangulation_tds.o $(LDFLAGS)

test_triangulation_3: test_triangulation_3.o 
	$(CGAL_CXX) -o test_triangulation_3 test_triangulation_3.o $(LDFLAGS)

test_delaunay_3: test_delaunay_3.o 
	$(CGAL_CXX) -o test_delaunay_3 test_delaunay_3.o $(LDFLAGS)

test_regular_3: test_regular_3.o 
	$(CGAL_CXX) -o test_regular_3 test_regular_3.o $(LDFLAGS)
clean:
	- rm -f *.o test_triangulation_3 test_triangulation_geom_traits \
	            test_triangulation_tds test_delaunay_3 \
		    test_regular_3 \
	            Program* Error* Compiler* error* Test*_triangulation_IO_3 \
	            Test_tds_IO_3 core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX)  $(CXXFLAGS) -c $<
