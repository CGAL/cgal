# This is the makefile for compiling a CGAL application.  
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	-I../../include \
	-I./include \
	 $(CGAL_CXXFLAGS) \
	 $(TESTSUITE_CXXFLAGS) \
	-W -Wall \
	-DCGAL_CHECK_EXPENSIVE 

	# -g
	# -DCGAL_NO_DEPRECATED_CODE \

	# -fprofile-arcs -ftest-coverage # pour Gcov
	# -Wold-style-cast
	# -pedantic
	# -DCGAL_CHECK_EXPENSIVE \
	# -DCGAL_CHECK_EXACTNESS \

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: test_triangulation_tds test_triangulation_3 test_delaunay_3 \
     test_regular_3 test_delaunay_hierarchy_3 test_delaunay_3_old_remove
	./test_triangulation_tds
	./test_triangulation_3
	./test_delaunay_3
	./test_delaunay_3_old_remove
	./test_delaunay_hierarchy_3
	./test_regular_3

test_triangulation_tds: test_triangulation_tds$(OBJ_EXT)  
	$(CGAL_CXX)  -o test_triangulation_tds test_triangulation_tds$(OBJ_EXT)  $(LDFLAGS)

test_triangulation_3: test_triangulation_3$(OBJ_EXT)  
	$(CGAL_CXX) -o test_triangulation_3 test_triangulation_3$(OBJ_EXT)  $(LDFLAGS)

test_delaunay_3: test_delaunay_3$(OBJ_EXT)  
	$(CGAL_CXX) -o test_delaunay_3 test_delaunay_3$(OBJ_EXT)  $(LDFLAGS)

test_delaunay_3_old_remove: test_delaunay_3_old_remove$(OBJ_EXT)  
	$(CGAL_CXX) -o test_delaunay_3_old_remove test_delaunay_3_old_remove$(OBJ_EXT)  $(LDFLAGS)

test_delaunay_hierarchy_3: test_delaunay_hierarchy_3$(OBJ_EXT)  
	$(CGAL_CXX) -o test_delaunay_hierarchy_3 test_delaunay_hierarchy_3$(OBJ_EXT)  $(LDFLAGS)

test_regular_3: test_regular_3$(OBJ_EXT)  
	$(CGAL_CXX) -o test_regular_3 test_regular_3$(OBJ_EXT)  $(LDFLAGS)

clean:
	- rm -f *$(OBJ_EXT)  test_triangulation_3 \
	            test_triangulation_tds test_delaunay_3 test_delaunay_3_old_remove \
		    test_delaunay_hierarchy_3 test_regular_3 \
	            Program* Error* Compiler* error* Test*_triangulation_IO_3 \
	            Test_tds_IO_3 core *.da *.gcov *.bb *.bbg

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT) :
	$(CGAL_CXX)  $(CXXFLAGS) -c $<
