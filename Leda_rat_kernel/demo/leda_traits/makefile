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
	   -I. -I../../include \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CGAL_GEOWIN_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#
#CGAL_poly_triang_2 \

all:            \
CGAL_alpha_shape_2 \
CGAL_arrangement_2 \
CGAL_constrained_triang_2 \
CGAL_convex_hull_2 \
CGAL_convex_hull_3 \
CGAL_delaunay_leda_2 \
CGAL_delaunay_leda_3 \
CGAL_delaunay_leda_event_3 \
CGAL_partition_leda_2 \
CGAL_segment_intersection_2 \
alpha_shape_2_leda \
alpha_shape_3_leda \
constrained_wi_2_leda \
convex_hull_2_leda \
delau_2_leda \
delau_3_kernel_check \
delau_3_leda \
delau_3_leda_counting \
delau_3_leda_new_check \
delau_3_leda_speed \
greene_2_leda \
nearest_neighbor \
optimal_convex_demo_leda \
partition_2_leda \
planar_map_example1 \
planar_map_example4 \
polygon_2_leda \
polyhedron_prog_incr_builder_leda \
polyhedron_prog_simple_leda \
quickhull_3_leda \
triang_2_leda \
event_test \
jarvis_anim_special \
graham_anim_special \
triangulation_anim

constrained_wi_2_leda$(EXE_EXT): constrained_wi_2_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained_wi_2_leda constrained_wi_2_leda$(OBJ_EXT) $(LDFLAGS)

triang_2_leda$(EXE_EXT): triang_2_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triang_2_leda triang_2_leda$(OBJ_EXT) $(LDFLAGS)
	
triang_2_leda$(EXE_EXT): triang_2_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triang_2_leda triang_2_leda$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_incr_builder_leda$(EXE_EXT): polyhedron_prog_incr_builder_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_incr_builder_leda polyhedron_prog_incr_builder_leda$(OBJ_EXT) $(LDFLAGS)

polyhedron_prog_simple_leda$(EXE_EXT): polyhedron_prog_simple_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_prog_simple_leda polyhedron_prog_simple_leda$(OBJ_EXT) $(LDFLAGS)

alpha_shape_2_leda$(EXE_EXT): alpha_shape_2_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)alpha_shape_2_leda alpha_shape_2_leda$(OBJ_EXT) $(LDFLAGS)

alpha_shape_3_leda$(EXE_EXT): alpha_shape_3_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)alpha_shape_3_leda alpha_shape_3_leda$(OBJ_EXT) $(LDFLAGS)

optimal_convex_demo_leda$(EXE_EXT): optimal_convex_demo_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)optimal_convex_demo_leda optimal_convex_demo_leda$(OBJ_EXT) $(LDFLAGS)

greene_2_leda$(EXE_EXT): greene_2_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)greene_2_leda greene_2_leda$(OBJ_EXT) $(LDFLAGS)

CGAL_partition_leda_2$(EXE_EXT): CGAL_partition_leda_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_partition_leda_2 CGAL_partition_leda_2$(OBJ_EXT) $(LDFLAGS)

CGAL_delaunay_leda_2$(EXE_EXT): CGAL_delaunay_leda_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_delaunay_leda_2 CGAL_delaunay_leda_2$(OBJ_EXT) $(LDFLAGS)

CGAL_delaunay_leda_3$(EXE_EXT): CGAL_delaunay_leda_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_delaunay_leda_3 CGAL_delaunay_leda_3$(OBJ_EXT) $(LDFLAGS)

CGAL_delaunay_leda_event_3$(EXE_EXT): CGAL_delaunay_leda_event_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_delaunay_leda_event_3 CGAL_delaunay_leda_event_3$(OBJ_EXT) $(LDFLAGS)

CGAL_poly_triang_2$(EXE_EXT): CGAL_poly_triang_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_poly_triang_2 CGAL_poly_triang_2$(OBJ_EXT) $(LDFLAGS)

CGAL_convex_hull_2$(EXE_EXT): CGAL_convex_hull_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_convex_hull_2 CGAL_convex_hull_2$(OBJ_EXT) $(LDFLAGS)

CGAL_convex_hull_3$(EXE_EXT): CGAL_convex_hull_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_convex_hull_3 CGAL_convex_hull_3$(OBJ_EXT) $(LDFLAGS)
	
CGAL_constrained_triang_2$(EXE_EXT): CGAL_constrained_triang_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_constrained_triang_2 CGAL_constrained_triang_2$(OBJ_EXT) $(LDFLAGS)

CGAL_alpha_shape_2$(EXE_EXT): CGAL_alpha_shape_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_alpha_shape_2 CGAL_alpha_shape_2$(OBJ_EXT) $(LDFLAGS)

CGAL_arrangment_2$(EXE_EXT): CGAL_arrangment_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_arrangement_2 CGAL_arrangement_2$(OBJ_EXT) $(LDFLAGS)

CGAL_segment_intersection_2$(EXE_EXT): CGAL_segment_intersection_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL_segment_intersection_2 CGAL_segment_intersection_2$(OBJ_EXT) $(LDFLAGS)

partition_2_leda$(EXE_EXT): partition_2_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)partition_2_leda partition_2_leda$(OBJ_EXT) $(LDFLAGS)

convex_hull_2_leda$(EXE_EXT): convex_hull_2_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)convex_hull_2_leda convex_hull_2_leda$(OBJ_EXT) $(LDFLAGS)

polygon_2_leda$(EXE_EXT): polygon_2_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polygon_2_leda polygon_2_leda$(OBJ_EXT) $(LDFLAGS)

nearest_neighbor$(EXE_EXT): nearest_neighbor$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nearest_neighbor nearest_neighbor$(OBJ_EXT) $(LDFLAGS)

delau_2_leda$(EXE_EXT): delau_2_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delau_2_leda delau_2_leda$(OBJ_EXT) $(LDFLAGS)

quickhull_3_leda$(EXE_EXT): quickhull_3_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)quickhull_3_leda quickhull_3_leda$(OBJ_EXT) $(LDFLAGS)

delau_3_leda$(EXE_EXT): delau_3_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delau_3_leda delau_3_leda$(OBJ_EXT) $(LDFLAGS)

delau_3_leda_counting$(EXE_EXT): delau_3_leda_counting$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delau_3_leda_counting delau_3_leda_counting$(OBJ_EXT) $(LDFLAGS)

delau_3_leda_new_check$(EXE_EXT): delau_3_leda_new_check$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delau_3_leda_new_check delau_3_leda_new_check$(OBJ_EXT) $(LDFLAGS)

delau_3_kernel_check$(EXE_EXT): delau_3_kernel_check$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delau_3_kernel_check delau_3_kernel_check$(OBJ_EXT) $(LDFLAGS)

delau_3_leda_speed$(EXE_EXT): delau_3_leda_speed$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delau_3_leda_speed delau_3_leda_speed$(OBJ_EXT) $(LDFLAGS)
	
triang_2_leda$(EXE_EXT): triang_2_leda$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triang_2_leda triang_2_leda$(OBJ_EXT) $(LDFLAGS)

planar_map_example1$(EXE_EXT): planar_map_example1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)planar_map_example1 planar_map_example1$(OBJ_EXT) $(LDFLAGS)

planar_map_example4$(EXE_EXT): planar_map_example4$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)planar_map_example4 planar_map_example4$(OBJ_EXT) $(LDFLAGS)

event_test$(EXE_EXT): event_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)event_test event_test$(OBJ_EXT) $(LDFLAGS)

jarvis_anim_special$(EXE_EXT): jarvis_anim_special$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)jarvis_anim_special jarvis_anim_special$(OBJ_EXT) $(LDFLAGS)

graham_anim_special$(EXE_EXT): graham_anim_special$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)graham_anim_special graham_anim_special$(OBJ_EXT) $(LDFLAGS)

triangulation_anim$(EXE_EXT): triangulation_anim$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangulation_anim triangulation_anim$(OBJ_EXT) $(LDFLAGS)

clean: \
CGAL_alpha_shape_2.clean \
CGAL_arrangement_2.clean \
CGAL_constrained_triang_2.clean \
CGAL_convex_hull_2.clean \
CGAL_convex_hull_3.clean \
CGAL_delaunay_leda_2.clean \
CGAL_delaunay_leda_3.clean \
CGAL_delaunay_leda_event_3.clean \
CGAL_partition_leda_2.clean \
CGAL_poly_triang_2.clean \
CGAL_segment_intersection_2.clean \
alpha_shape_2_leda.clean \
alpha_shape_3_leda.clean \
constrained_wi_2_leda.clean \
convex_hull_2_leda.clean \
delau_2_leda.clean \
delau_3_kernel_check.clean \
delau_3_leda.clean \
delau_3_leda_counting.clean \
delau_3_leda_new_check.clean \
delau_3_leda_speed.clean \
greene_2_leda.clean \
nearest_neighbor.clean \
optimal_convex_demo_leda.clean \
partition_2_leda.clean \
planar_map_example1.clean \
planar_map_example4.clean \
polygon_2_leda.clean \
polyhedron_prog_incr_builder_leda.clean \
polyhedron_prog_simple_leda.clean \
quickhull_3_leda.clean \
triang_2_leda.clean \
event_test.clean \
graham_anim_special.clean \
triangulation_anim.clean \
CGAL_delaunay_leda_event_3.clean \
jarvis_anim_special.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

