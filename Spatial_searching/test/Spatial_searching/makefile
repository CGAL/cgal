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
           -I../../include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

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
                Building_kd_tree_with_ddim_points$(EXE_EXT) \
                Building_kd_tree_with_own_pointtype$(EXE_EXT) \
                Circular_query$(EXE_EXT) \
                Iso_rectangle_2_query$(EXE_EXT) \
                K_neighbor_search_manhattan_distance_isobox_point$(EXE_EXT) \
                K_neighbor_search_with_circle$(EXE_EXT) \
                Orthogonal_incremental_neighbor_search$(EXE_EXT) \
                Orthogonal_k_neighbor_search$(EXE_EXT) \
                Range_searching$(EXE_EXT) \
                Splitters$(EXE_EXT) 

Building_kd_tree_with_ddim_points$(EXE_EXT): Building_kd_tree_with_ddim_points$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Building_kd_tree_with_ddim_points Building_kd_tree_with_ddim_points$(OBJ_EXT) $(LDFLAGS)

Building_kd_tree_with_own_pointtype$(EXE_EXT): Building_kd_tree_with_own_pointtype$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Building_kd_tree_with_own_pointtype Building_kd_tree_with_own_pointtype$(OBJ_EXT) $(LDFLAGS)

Circular_query$(EXE_EXT): Circular_query$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Circular_query Circular_query$(OBJ_EXT) $(LDFLAGS)

Iso_rectangle_2_query$(EXE_EXT): Iso_rectangle_2_query$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Iso_rectangle_2_query Iso_rectangle_2_query$(OBJ_EXT) $(LDFLAGS)

K_neighbor_search_manhattan_distance_isobox_point$(EXE_EXT): K_neighbor_search_manhattan_distance_isobox_point$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)K_neighbor_search_manhattan_distance_isobox_point K_neighbor_search_manhattan_distance_isobox_point$(OBJ_EXT) $(LDFLAGS)

K_neighbor_search_with_circle$(EXE_EXT): K_neighbor_search_with_circle$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)K_neighbor_search_with_circle K_neighbor_search_with_circle$(OBJ_EXT) $(LDFLAGS)

Orthogonal_incremental_neighbor_search$(EXE_EXT): Orthogonal_incremental_neighbor_search$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Orthogonal_incremental_neighbor_search Orthogonal_incremental_neighbor_search$(OBJ_EXT) $(LDFLAGS)

Orthogonal_k_neighbor_search$(EXE_EXT): Orthogonal_k_neighbor_search$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Orthogonal_k_neighbor_search Orthogonal_k_neighbor_search$(OBJ_EXT) $(LDFLAGS)

Range_searching$(EXE_EXT): Range_searching$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Range_searching Range_searching$(OBJ_EXT) $(LDFLAGS)

Splitters$(EXE_EXT): Splitters$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Splitters Splitters$(OBJ_EXT) $(LDFLAGS)

clean: \
                   Building_kd_tree_with_ddim_points.clean \
                   Building_kd_tree_with_own_pointtype.clean \
                   Circular_query.clean \
                   Iso_rectangle_2_query.clean \
                   K_neighbor_search_manhattan_distance_isobox_point.clean \
                   K_neighbor_search_with_circle.clean \
                   Orthogonal_incremental_neighbor_search.clean \
                   Orthogonal_k_neighbor_search.clean \
                   Range_searching.clean \
                   Splitters.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

