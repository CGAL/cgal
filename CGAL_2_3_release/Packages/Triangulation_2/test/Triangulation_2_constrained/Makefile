# This is the makefile for compiling a CGAL application.  
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.


CGAL_MAKEFILE =/0/prisme_util/CGAL/CGAL-I/make/makefile_i586_Linux-2.2.14_g++-2.95.3_LEDA
include $(CGAL_MAKEFILE)


#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	-I../Triangulation_2/include \
	-I./include \
	-I../../include \
	$(TESTSUITE_CXXFLAGS) \
	$(CGAL_CXXFLAGS)\
	-g 


##To test the Kernel
#	-I../../../Kernel/C2/include \
#	-I../../../Kernel/C3/include \
#	-I../../../Kernel/Cartesian_basic/include \





#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)
#	  $(LONG_NAME_PROBLEM_LDFLAGS) \

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#


test_constrained_triangulation_2: test_constrained_triangulation_2.o	
	$(CGAL_CXX) -o test_constrained_triangulation_2 test_constrained_triangulation_2.o $(LDFLAGS)

test_const_del_triangulation_2: test_const_del_triangulation_2.o
	$(CGAL_CXX) -o test_const_del_triangulation_2 test_const_del_triangulation_2.o $(LDFLAGS)

test_delaunay_hierarchy_2: test_delaunay_hierarchy_2.o
	$(CGAL_CXX) -o test_delaunay_hierarchy_2 test_delaunay_hierarchy_2.o $(LDFLAGS)

clean:
	- rm -f *.o  *.obj *.exe \
	test_constrained_triangulation_2 \
	test_const_del_triangulation_2 \
	test_delaunay_hierarchy_2 \
	 core \
	T??.triangulation  file_tds? 

realclean: clean
	- rm -f Program* Error* Compiler* error*

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX)  $(CXXFLAGS) -c $<
