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
	-I../Triangulation_2/include \
	$(TESTSUITE_CXXFLAGS) \
	$(EXTRA_FLAGS) \
	$(CGAL_CXXFLAGS) \
	$(LONG_NAME_PROBLEM_CXXFLAGS) \
	$(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
	$(TESTSUITE_LIBPATH) \
	$(CGAL_LIBPATH)

LDFLAGS = \
	$(TESTSUITE_LDFLAGS) \
	$(LONG_NAME_PROBLEM_LDFLAGS) \
	$(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: \
		test_delaunay_hierarchy_2$(EXE_EXT) \
		test_delaunay_triangulation_2$(EXE_EXT) \
		test_regular_triangulation_2$(EXE_EXT) 


test_delaunay_hierarchy_2$(EXE_EXT): test_delaunay_hierarchy_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_delaunay_hierarchy_2 test_delaunay_hierarchy_2$(OBJ_EXT) $(LDFLAGS)

test_delaunay_triangulation_2$(EXE_EXT): test_delaunay_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_delaunay_triangulation_2 test_delaunay_triangulation_2$(OBJ_EXT) $(LDFLAGS)

test_regular_triangulation_2$(EXE_EXT): test_regular_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_regular_triangulation_2 test_regular_triangulation_2$(OBJ_EXT) $(LDFLAGS)

clean: \
		test_delaunay_hierarchy_2.clean \
		test_delaunay_triangulation_2.clean \
		test_regular_triangulation_2.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

