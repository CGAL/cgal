# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

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

all:            \
                all_furthest_neighbors_2_test$(EXE_EXT) \
                extremal_polygon_2_test$(EXE_EXT) \
                rectangular_p_center_2_random1_test$(EXE_EXT) \
                rectangular_p_center_2_random2_test$(EXE_EXT) \
                sorted_matrix_search_test$(EXE_EXT) 

all_furthest_neighbors_2_test$(EXE_EXT): all_furthest_neighbors_2_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)all_furthest_neighbors_2_test all_furthest_neighbors_2_test$(OBJ_EXT) $(LDFLAGS)

extremal_polygon_2_test$(EXE_EXT): extremal_polygon_2_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)extremal_polygon_2_test extremal_polygon_2_test$(OBJ_EXT) $(LDFLAGS)

rectangular_p_center_2_random1_test$(EXE_EXT): rectangular_p_center_2_random1_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rectangular_p_center_2_random1_test rectangular_p_center_2_random1_test$(OBJ_EXT) $(LDFLAGS)

rectangular_p_center_2_random2_test$(EXE_EXT): rectangular_p_center_2_random2_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rectangular_p_center_2_random2_test rectangular_p_center_2_random2_test$(OBJ_EXT) $(LDFLAGS)

sorted_matrix_search_test$(EXE_EXT): sorted_matrix_search_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)sorted_matrix_search_test sorted_matrix_search_test$(OBJ_EXT) $(LDFLAGS)

clean: \
                   all_furthest_neighbors_2_test.clean \
                   extremal_polygon_2_test.clean \
                   rectangular_p_center_2_random1_test.clean \
                   rectangular_p_center_2_random2_test.clean \
                   sorted_matrix_search_test.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

