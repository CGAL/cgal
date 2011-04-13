# This is the makefile for compiling a CGAL application.  
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# pour purify, virer LEDA
#CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-I/make/makefile_sparc_SunOS-5.6_g++-2.95.3_LEDA
#CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-I/make/makefile_i686_Linux-2.2.7_g++-2.95.1_LEDA
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

# pour purify :
#CGAL_CXX=purify -best-effort /usr/local/gcc/bin/g++
# pour insure :
# CGAL_CXX=insure

CXXFLAGS = \
	-I../../include \
	-I./include \
	 $(CGAL_CXXFLAGS) \
	 $(TESTSUITE_CXXFLAGS) \
	-W -Wall \
	-DCGAL_CHECK_EXPENSIVE \
	-g # -Wold-style-cast # -pedantic
        # -I$(HOME)/prog/Cgal/Kernel_basic/include \
        # -I$(HOME)/prog/Cgal/Cartesian_kernel/include \
	# -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS \
	# -DCGAL_NO_ASSERTIONS -DCGAL_NO_WARNINGS \
	# -DCGAL_CHECK_EXPENSIVE \
	# -B$(UTIL)/Binutils/$(UNAME)/bin/ \
	# -gdwarf-2
	# -O2 -fno-default-inline \

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

# pour purify :
#CGAL_LDFLAGS= $(CGAL_LIB_DIR)/$(CGAL_OS_COMPILER)/libCGAL.a \
              # -L$(GMP_LIB_DIR) -lgmp -lm 

LDFLAGS = $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)
#          -B$(UTIL)/Binutils/$(UNAME)/bin/ \

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: test_triangulation_tds \
     test_triangulation_3 test_delaunay_3 \
     test_regular_3 test_delaunay_hierarchy_3
	./test_triangulation_tds
	./test_triangulation_3
	./test_delaunay_3
	./test_delaunay_hierarchy_3
	./test_regular_3 
#        -B$(UTIL)/Binutils/$(UNAME)/bin/ \

test_triangulation_tds: ../../src/Triangulation_3.o test_triangulation_tds.o 
	$(CGAL_CXX)  -o test_triangulation_tds Triangulation_3.o test_triangulation_tds.o $(LDFLAGS)

test_triangulation_3: ../../src/Triangulation_3.o  test_triangulation_3.o 
	$(CGAL_CXX) -o test_triangulation_3 Triangulation_3.o test_triangulation_3.o $(LDFLAGS)

test_delaunay_3: ../../src/Triangulation_3.o  test_delaunay_3.o 
	$(CGAL_CXX) -o test_delaunay_3 Triangulation_3.o test_delaunay_3.o $(LDFLAGS)

test_delaunay_hierarchy_3: ../../src/Triangulation_3.o test_delaunay_hierarchy_3.o 
	$(CGAL_CXX) -o test_delaunay_hierarchy_3 Triangulation_3.o test_delaunay_hierarchy_3.o $(LDFLAGS)

test_regular_3: ../../src/Triangulation_3.o  test_regular_3.o 
	$(CGAL_CXX) -o test_regular_3 Triangulation_3.o test_regular_3.o $(LDFLAGS)

clean:
	- rm -f *.o test_triangulation_3 \
	            test_triangulation_tds test_delaunay_3 \
		    test_delaunay_hierarchy_3 test_regular_3 \
	            Program* Error* Compiler* error* Test*_triangulation_IO_3 \
	            Test_tds_IO_3 core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX)  $(CXXFLAGS) -c $<
