
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = 

CGDL = \
          $(LDFLAGS) $(CGAL_GEOWIN_LDFLAGS)

CGAL_alpha_shape_2$(EXE_EXT) :  CGAL_alpha_shape_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_alpha_shape_2 CGAL_alpha_shape_2$(OBJ_EXT)  $(CGDL)

CGAL_arrangement_2$(EXE_EXT) :  CGAL_arrangement_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_arrangement_2 CGAL_arrangement_2$(OBJ_EXT)  $(CGDL)

CGAL_geometry_2$(EXE_EXT) :  CGAL_geometry_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_geometry_2 CGAL_geometry_2$(OBJ_EXT)  $(CGDL)
	
CGAL_delaunay_triang_3$(EXE_EXT) :  CGAL_delaunay_triang_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_delaunay_triang_3 CGAL_delaunay_triang_3$(OBJ_EXT)  $(CGDL)

CGAL_convex_hull_2$(EXE_EXT) :  CGAL_convex_hull_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_convex_hull_2 CGAL_convex_hull_2$(OBJ_EXT)  $(CGDL)

CGAL_kdtree_2$(EXE_EXT) :  CGAL_kdtree_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_kdtree_2 CGAL_kdtree_2$(OBJ_EXT)  $(CGDL)

CGAL_kdtree_3$(EXE_EXT) :  CGAL_kdtree_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_kdtree_3 CGAL_kdtree_3$(OBJ_EXT)  $(CGDL)

CGAL_convex_hull_3$(EXE_EXT) :  CGAL_convex_hull_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_convex_hull_3 CGAL_convex_hull_3$(OBJ_EXT)  $(CGDL)

CGAL_segment_intersection_2$(EXE_EXT) :  CGAL_segment_intersection_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_segment_intersection_2 CGAL_segment_intersection_2$(OBJ_EXT)  $(CGDL)

CGAL_polygon_2$(EXE_EXT) :  CGAL_polygon_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_polygon_2 CGAL_polygon_2$(OBJ_EXT)  $(CGDL)

CGAL_max_k_gon_2$(EXE_EXT) :  CGAL_max_k_gon_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_max_k_gon_2 CGAL_max_k_gon_2$(OBJ_EXT)  $(CGDL)

CGAL_all_furthest_neighbors_2$(EXE_EXT) :  CGAL_all_furthest_neighbors_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_all_furthest_neighbors_2 CGAL_all_furthest_neighbors_2$(OBJ_EXT)  $(CGDL)

CGAL_demo$(EXE_EXT) :  CGAL_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_demo CGAL_demo$(OBJ_EXT)  $(CGDL)

CGAL_regular_triang_2$(EXE_EXT) :  CGAL_regular_triang_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_regular_triang_2 CGAL_regular_triang_2$(OBJ_EXT)  $(CGDL)
	
CGAL_triangulation_2$(EXE_EXT) :  CGAL_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_triangulation_2 CGAL_triangulation_2$(OBJ_EXT)  $(CGDL)

CGAL_constrained_triang_2$(EXE_EXT) :  CGAL_constrained_triang_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_constrained_triang_2 CGAL_constrained_triang_2$(OBJ_EXT)  $(CGDL)
		
CGAL_delaunay_triang_2$(EXE_EXT) :  CGAL_delaunay_triang_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_delaunay_triang_2 CGAL_delaunay_triang_2$(OBJ_EXT)  $(CGDL)
	
CGAL_min_circle_2$(EXE_EXT) :  CGAL_min_circle_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_min_circle_2 CGAL_min_circle_2$(OBJ_EXT)  $(CGDL)

CGAL_min_ellipse_2$(EXE_EXT) :  CGAL_min_ellipse_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_min_ellipse_2 CGAL_min_ellipse_2$(OBJ_EXT)  $(CGDL)

CGAL_min_sphere_3$(EXE_EXT) :  CGAL_min_sphere_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_min_sphere_3 CGAL_min_sphere_3$(OBJ_EXT)  $(CGDL)

CGAL_planar_map_2$(EXE_EXT) :  CGAL_planar_map_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_planar_map_2 CGAL_planar_map_2$(OBJ_EXT)  $(CGDL)

all: CGAL_alpha_shape_2 CGAL_arrangement_2 CGAL_constrained_triang_2 CGAL_min_sphere_3 CGAL_regular_triang_2 CGAL_delaunay_triang_3 CGAL_geometry_2 CGAL_min_circle_2 CGAL_delaunay_triang_2 CGAL_demo CGAL_convex_hull_2 CGAL_segment_intersection_2 CGAL_triangulation_2 CGAL_max_k_gon_2 CGAL_all_furthest_neighbors_2 CGAL_polygon_2 CGAL_convex_hull_3 CGAL_kdtree_2 CGAL_kdtree_3 CGAL_min_ellipse_2


clean : \
	CGAL_alpha_shape_2.clean CGAL_arrangement_2.clean \
	 CGAL_constrained_triang_2.clean CGAL_min_sphere_3.clean \
	 CGAL_regular_triang_2.clean CGAL_delaunay_triang_3.clean \
	 CGAL_geometry_2.clean CGAL_min_circle_2.clean \
	 CGAL_delaunay_triang_2.clean CGAL_demo.clean \ 
	CGAL_convex_hull_2.clean  CGAL_segment_intersection_2.clean \
	 CGAL_triangulation_2.clean CGAL_max_k_gon_2.clean \ 
	 CGAL_all_furthest_neighbors_2.clean CGAL_polygon_2.clean \ 
	 CGAL_convex_hull_3.clean CGAL_kdtree_2.clean \ 
	 CGAL_kdtree_3.clean CGAL_min_ellipse_2.clean \ 
	 CGAL_planar_map_2.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<




