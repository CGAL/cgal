# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = c:/cgal/CGAL-3.2-I-384/make/makefile_i686_CYGWINNT-5.1_CL.EXE-1310
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include \
           -I$(HOME)/CGAL/trunk/Kernel_23/include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	    -g

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

all:            \
                test_delaunay_3$(EXE_EXT) \
                test_delaunay_3_old_remove$(EXE_EXT) \
                test_delaunay_hierarchy_3$(EXE_EXT) \
                test_regular_3$(EXE_EXT) \
                test_regular_as_delaunay_3$(EXE_EXT) \
                test_regular_traits_3$(EXE_EXT) \
                test_regular_remove_3$(EXE_EXT) \
                test_triangulation_3$(EXE_EXT) \
                test_triangulation_tds$(EXE_EXT) 

test_delaunay_3$(EXE_EXT): test_delaunay_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_delaunay_3 test_delaunay_3$(OBJ_EXT) $(LDFLAGS)

test_delaunay_3_old_remove$(EXE_EXT): test_delaunay_3_old_remove$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_delaunay_3_old_remove test_delaunay_3_old_remove$(OBJ_EXT) $(LDFLAGS)

test_delaunay_hierarchy_3$(EXE_EXT): test_delaunay_hierarchy_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_delaunay_hierarchy_3 test_delaunay_hierarchy_3$(OBJ_EXT) $(LDFLAGS)

test_regular_3$(EXE_EXT): test_regular_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_regular_3 test_regular_3$(OBJ_EXT) $(LDFLAGS)

test_regular_as_delaunay_3$(EXE_EXT): test_regular_as_delaunay_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_regular_as_delaunay_3 test_regular_as_delaunay_3$(OBJ_EXT) $(LDFLAGS)

test_regular_traits_3$(EXE_EXT): test_regular_traits_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_regular_traits_3 test_regular_traits_3$(OBJ_EXT) $(LDFLAGS)

test_regular_remove_3$(EXE_EXT): test_regular_remove_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_regular_remove_3 test_regular_remove_3$(OBJ_EXT) $(LDFLAGS)

test_triangulation_3$(EXE_EXT): test_triangulation_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_triangulation_3 test_triangulation_3$(OBJ_EXT) $(LDFLAGS)

test_triangulation_tds$(EXE_EXT): test_triangulation_tds$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_triangulation_tds test_triangulation_tds$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_delaunay_3.clean \
                   test_delaunay_3_old_remove.clean \
                   test_delaunay_hierarchy_3.clean \
                   test_regular_3.clean \
                   test_regular_as_delaunay_3.clean \
                   test_regular_traits_3.clean \
                   test_regular_remove_3.clean \
                   test_triangulation_3.clean \
                   test_triangulation_tds.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

