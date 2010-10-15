# Created by the script cgal_create_makefile
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
           $(ADDITIONAL_CXXFLAGS) \
	   $(OPTIMIZATION_FLAGS) \
           $(KERNEL_FLAGS) \
           -DCGAL_NEF3_DUMP_STATISTICS \
           -DCGAL_NEF3_TIMER_INITIALIZE_KDTREE \
           -DCGAL_NEF3_TIMER_PLUECKER \
	   -DCGAL_NEF3_TIMER_PLANE_SWEEPS \
           -DCGAL_NEF3_TIMER_SPHERE_SWEEPS \
           -DCGAL_NEF3_TIMER_POINT_LOCATION \
           -DCGAL_NEF3_TIMER_SIMPLIFICATION \
           -DCGAL_NEF3_TIMER_INTERSECTION \
           -DCGAL_NEF3_TIMER_EXTERNAL_STRUCTURE \
           -DCGAL_NEF3_TIMER_OVERLAY \
	   -DCGAL_NEF3_TIMER_BINARY_OPERATION \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                union_of_rotated_cylinders$(EXE_EXT) \
                union_of_rotated_cylinders2$(EXE_EXT) \
	        two_tetrahedra_fields$(EXE_EXT) \
	       	tetrahedra_and_grid$(EXE_EXT) \
                big_and_small$(EXE_EXT) \
                quadratic$(EXE_EXT) \
                quadratic2$(EXE_EXT) \
                create_quadratic2$(EXE_EXT) \
                convert_to_sat$(EXE_EXT) \
                tetrahedra_as_sat$(EXE_EXT) \
                grid_as_sat$(EXE_EXT) \
                abdruck$(EXE_EXT) \
                complex_sphere_map$(EXE_EXT) \
                complex_facet$(EXE_EXT) \
                point_location$(EXE_EXT) \
		point_location2$(EXE_EXT) \
		create_complex_sphere_map_data$(EXE_EXT) \
		gaussian_map_overlay$(EXE_EXT) \
		move_x_half$(EXE_EXT) \
	  	binary_operation$(EXE_EXT) \
	  	minkowski_sum$(EXE_EXT) \
	  	incremental$(EXE_EXT)

union_of_rotated_cylinders$(EXE_EXT): union_of_rotated_cylinders$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)union_of_rotated_cylinders union_of_rotated_cylinders$(OBJ_EXT) $(LDFLAGS)
union_of_rotated_cylinders2$(EXE_EXT): union_of_rotated_cylinders2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)union_of_rotated_cylinders2 union_of_rotated_cylinders2$(OBJ_EXT) $(LDFLAGS)
two_tetrahedra_fields$(EXE_EXT): two_tetrahedra_fields$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)two_tetrahedra_fields two_tetrahedra_fields$(OBJ_EXT) $(LDFLAGS)
tetrahedra_and_grid$(EXE_EXT): tetrahedra_and_grid$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tetrahedra_and_grid tetrahedra_and_grid$(OBJ_EXT) $(LDFLAGS)
big_and_small$(EXE_EXT): big_and_small$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)big_and_small big_and_small$(OBJ_EXT) $(LDFLAGS)
quadratic$(EXE_EXT): quadratic$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)quadratic quadratic$(OBJ_EXT) $(LDFLAGS)
quadratic2$(EXE_EXT): quadratic2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)quadratic2 quadratic2$(OBJ_EXT) $(LDFLAGS)
create_quadratic2$(EXE_EXT): create_quadratic2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)create_quadratic2 create_quadratic2$(OBJ_EXT) $(LDFLAGS)
convert_to_sat$(EXE_EXT): convert_to_sat$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)convert_to_sat convert_to_sat$(OBJ_EXT) $(LDFLAGS)
tetrahedra_as_sat$(EXE_EXT): tetrahedra_as_sat$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tetrahedra_as_sat tetrahedra_as_sat$(OBJ_EXT) $(LDFLAGS)
grid_as_sat$(EXE_EXT): grid_as_sat$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)grid_as_sat grid_as_sat$(OBJ_EXT) $(LDFLAGS)
abdruck$(EXE_EXT): abdruck$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)abdruck abdruck$(OBJ_EXT) $(LDFLAGS)
complex_sphere_map$(EXE_EXT): complex_sphere_map$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)complex_sphere_map complex_sphere_map$(OBJ_EXT) $(LDFLAGS)
complex_facet$(EXE_EXT): complex_facet$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)complex_facet complex_facet$(OBJ_EXT) $(LDFLAGS)
point_location$(EXE_EXT): point_location$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)point_location point_location$(OBJ_EXT) $(LDFLAGS)
point_location2$(EXE_EXT): point_location2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)point_location2 point_location2$(OBJ_EXT) $(LDFLAGS)
create_complex_sphere_map_data$(EXE_EXT): create_complex_sphere_map_data$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)create_complex_sphere_map_data create_complex_sphere_map_data$(OBJ_EXT) $(LDFLAGS)
gaussian_map_overlay$(EXE_EXT): gaussian_map_overlay$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)gaussian_map_overlay gaussian_map_overlay$(OBJ_EXT) $(LDFLAGS)
move_x_half$(EXE_EXT): move_x_half$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)move_x_half move_x_half$(OBJ_EXT) $(LDFLAGS)
binary_operation$(EXE_EXT): binary_operation$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)binary_operation binary_operation$(OBJ_EXT) $(LDFLAGS)
minkowski_sum$(EXE_EXT): minkowski_sum$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)minkowski_sum minkowski_sum$(OBJ_EXT) $(LDFLAGS)
incremental$(EXE_EXT): incremental$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)incremental incremental$(OBJ_EXT) $(LDFLAGS)

clean: \
                   union_of_rotated_cylinders.clean \
                   union_of_rotated_cylinders2.clean \
	           two_tetrahedra_fields.clean \
	           tetrahedra_and_grid.clean \
	           big_and_small.clean \
	           quadratic.clean \
	           quadratic2.clean \
	           create_quadratic2.clean \
	           convert_to_sat.clean \
	           tetrahedra_as_sat.clean \
	           grid_as_sat.clean \
	           abdruck.clean \
	           complex_sphere_map.clean \
	           complex_facet.clean \
	           point_location.clean \
	           point_location2.clean \
	           create_complex_sphere_map_data.clean \
	           gaussian_map_overlay.clean \
                   move_x_half.clean \
	           binary_operation.clean \
	           minkowski_sum.clean \
	           incremental.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

