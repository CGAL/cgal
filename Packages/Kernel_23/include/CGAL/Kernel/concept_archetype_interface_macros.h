// ======================================================================
//
// Copyright (c) 1999,2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : concept_archetype_interface_macros.h
// package       : Kernel_23
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Matthias Baesken
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

// for the original version see Kernel/interface_macros.h

// This file is intentionally not protected against re-inclusion.
// It's aimed at being included from within a kernel traits class, this
// way we share more code.

// It is the responsability of the including file to correctly set the 2
// macros CGAL_Kernel_pred and CGAL_Kernel_cons.
// And they are #undefed at the end of this file.

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ANGLE_2)
CGAL_Kernel_pred(Angle_2,
		 angle_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ANGLE_3)
CGAL_Kernel_pred(Angle_3,
		 angle_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_ARE_ORDERED_ALONG_LINE_2)
CGAL_Kernel_pred(Are_ordered_along_line_2,
		 are_ordered_along_line_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_ARE_ORDERED_ALONG_LINE_3)
CGAL_Kernel_pred(Are_ordered_along_line_3,
		 are_ordered_along_line_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_ARE_STRICTLY_ORDERED_ALONG_LINE_2)
CGAL_Kernel_pred(Are_strictly_ordered_along_line_2,
		 are_strictly_ordered_along_line_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_ARE_STRICTLY_ORDERED_ALONG_LINE_3)
CGAL_Kernel_pred(Are_strictly_ordered_along_line_3,
		 are_strictly_ordered_along_line_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ASSIGN_2)
CGAL_Kernel_cons(Assign_2,
		 assign_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ASSIGN_3)
CGAL_Kernel_cons(Assign_3,
		 assign_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_BOUNDED_SIDE_2)
CGAL_Kernel_pred(Bounded_side_2,
		 bounded_side_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_BOUNDED_SIDE_3)
CGAL_Kernel_pred(Bounded_side_3,
		 bounded_side_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COLLINEAR_ARE_ORDERED_ALONG_LINE_2)
CGAL_Kernel_pred(Collinear_are_ordered_along_line_2,
		 collinear_are_ordered_along_line_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COLLINEAR_ARE_ORDERED_ALONG_LINE_3)
CGAL_Kernel_pred(Collinear_are_ordered_along_line_3,
		 collinear_are_ordered_along_line_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COLLINEAR_ARE_STRICTLY_ORDERED_ALONG_LINE_2)
CGAL_Kernel_pred(Collinear_are_strictly_ordered_along_line_2,
		 collinear_are_strictly_ordered_along_line_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COLLINEAR_ARE_STRICTLY_ORDERED_ALONG_LINE_3)
CGAL_Kernel_pred(Collinear_are_strictly_ordered_along_line_3,
		 collinear_are_strictly_ordered_along_line_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COLLINEAR_HAS_ON_2)
CGAL_Kernel_pred(Collinear_has_on_2,
		 collinear_has_on_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COLLINEAR_2)
CGAL_Kernel_pred(Collinear_2,
		 collinear_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COLLINEAR_3)
CGAL_Kernel_pred(Collinear_3,
		 collinear_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPARE_ANGLE_WITH_X_AXIS_2)
CGAL_Kernel_pred(Compare_angle_with_x_axis_2,
		 compare_angle_with_x_axis_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPARE_DISTANCE_2)
CGAL_Kernel_pred(Compare_distance_2,
		 compare_distance_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPARE_DISTANCE_3)
CGAL_Kernel_pred(Compare_distance_3,
		 compare_distance_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPARE_SLOPE_2)
CGAL_Kernel_pred(Compare_slope_2,
		 compare_slope_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPARE_X_AT_Y_2)
CGAL_Kernel_pred(Compare_x_at_y_2,
		 compare_x_at_y_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPARE_XYZ_3)
CGAL_Kernel_pred(Compare_xyz_3,
		 compare_xyz_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPARE_XY_2)
CGAL_Kernel_pred(Compare_xy_2,
		 compare_xy_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPARE_XY_3)
CGAL_Kernel_pred(Compare_xy_3,
		 compare_xy_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPARE_X_2)
CGAL_Kernel_pred(Compare_x_2,
		 compare_x_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPARE_X_3)
CGAL_Kernel_pred(Compare_x_3,
		 compare_x_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPARE_Y_AT_X_2)
CGAL_Kernel_pred(Compare_y_at_x_2,
		 compare_y_at_x_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPARE_Y_2)
CGAL_Kernel_pred(Compare_y_2,
		 compare_y_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPARE_Y_3)
CGAL_Kernel_pred(Compare_y_3,
		 compare_y_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPARE_Z_3)
CGAL_Kernel_pred(Compare_z_3,
		 compare_z_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPUTE_AREA_2)
CGAL_Kernel_cons(Compute_area_2,
		 compute_area_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COMPUTE_AREA_3)
CGAL_Kernel_cons(Compute_squared_area_3,
		 compute_squared_area_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPUTE_SQUARED_DISTANCE_2)
CGAL_Kernel_cons(Compute_squared_distance_2,
		 compute_squared_distance_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPUTE_SQUARED_DISTANCE_3)
CGAL_Kernel_cons(Compute_squared_distance_3,
		 compute_squared_distance_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPUTE_SQUARED_LENGTH_2)
CGAL_Kernel_cons(Compute_squared_length_2,
		 compute_squared_length_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPUTE_SQUARED_LENGTH_3)
CGAL_Kernel_cons(Compute_squared_length_3,
		 compute_squared_length_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPUTE_SQUARED_RADIUS_2)
CGAL_Kernel_cons(Compute_squared_radius_2,
		 compute_squared_radius_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPUTE_SQUARED_RADIUS_3)
CGAL_Kernel_cons(Compute_squared_radius_3,
		 compute_squared_radius_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COMPUTE_VOLUME_3)
CGAL_Kernel_cons(Compute_volume_3,
		 compute_volume_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_BASE_VECTOR_3)
CGAL_Kernel_cons(Construct_base_vector_3,
		 construct_base_vector_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_BISECTOR_2)
CGAL_Kernel_cons(Construct_bisector_2,
		 construct_bisector_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_CENTER_2)
CGAL_Kernel_cons(Construct_center_2,
		 construct_center_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_CENTER_3)
CGAL_Kernel_cons(Construct_center_3,
		 construct_center_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_CENTROID_2)
CGAL_Kernel_cons(Construct_centroid_2,
		 construct_centroid_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_CENTROID_3)
CGAL_Kernel_cons(Construct_centroid_3,
		 construct_centroid_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_CIRCLE_2)
CGAL_Kernel_cons(Construct_circle_2,
		 construct_circle_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_CIRCUMCENTER_2)
CGAL_Kernel_cons(Construct_circumcenter_2,
		 construct_circumcenter_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_CIRCUMCENTER_3)
CGAL_Kernel_cons(Construct_circumcenter_3,
		 construct_circumcenter_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_CROSS_PRODUCT_VECTOR_3)
CGAL_Kernel_cons(Construct_cross_product_vector_3,
		 construct_cross_product_vector_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_DIRECTION_2)
CGAL_Kernel_cons(Construct_direction_2,
		 construct_direction_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_DIRECTION_3)
CGAL_Kernel_cons(Construct_direction_3,
		 construct_direction_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_ISO_CUBOID_3)
CGAL_Kernel_cons(Construct_iso_cuboid_3,
		 construct_iso_cuboid_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_ISO_RECTANGLE_2)
CGAL_Kernel_cons(Construct_iso_rectangle_2,
		 construct_iso_rectangle_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_LIFTED_POINT_3)
CGAL_Kernel_cons(Construct_lifted_point_3,
		 construct_lifted_point_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CONSTRUCT_LINE_2)
CGAL_Kernel_cons(Construct_line_2,
		 construct_line_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CONSTRUCT_LINE_3)
CGAL_Kernel_cons(Construct_line_3,
		 construct_line_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_MIDPOINT_2)
CGAL_Kernel_cons(Construct_midpoint_2,
		 construct_midpoint_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_MIDPOINT_3)
CGAL_Kernel_cons(Construct_midpoint_3,
		 construct_midpoint_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OBJECT_2)
CGAL_Kernel_cons(Construct_object_2,
                 construct_object_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CONSTRUCT_OBJECT_3)
CGAL_Kernel_cons(Construct_object_3,
                 construct_object_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_CIRCLE_2)
CGAL_Kernel_cons(Construct_opposite_circle_2,
		 construct_opposite_circle_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_DIRECTION_2)
CGAL_Kernel_cons(Construct_opposite_direction_2,
		 construct_opposite_direction_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_DIRECION_3)
CGAL_Kernel_cons(Construct_opposite_direction_3,
		 construct_opposite_direction_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_LINE_2)
CGAL_Kernel_cons(Construct_opposite_line_2,
		 construct_opposite_line_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_LINE_3)
CGAL_Kernel_cons(Construct_opposite_line_3,
		 construct_opposite_line_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_PLANE_3)
CGAL_Kernel_cons(Construct_opposite_plane_3,
		 construct_opposite_plane_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_RAY_2)
CGAL_Kernel_cons(Construct_opposite_ray_2,
		 construct_opposite_ray_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_RAY_3)
CGAL_Kernel_cons(Construct_opposite_ray_3,
		 construct_opposite_ray_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_SEGMENT_2)
CGAL_Kernel_cons(Construct_opposite_segment_2,
		 construct_opposite_segment_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_SEGMENT_3)
CGAL_Kernel_cons(Construct_opposite_segment_3,
		 construct_opposite_segment_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_SPHERE_3)
CGAL_Kernel_cons(Construct_opposite_sphere_3,
		 construct_opposite_sphere_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_TRIANGLE_2)
CGAL_Kernel_cons(Construct_opposite_triangle_2,
		 construct_opposite_triangle_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_VECTOR_2)
CGAL_Kernel_cons(Construct_opposite_vector_2,
		 construct_opposite_vector_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_OPPOSITE_VECTOR_3)
CGAL_Kernel_cons(Construct_opposite_vector_3,
		 construct_opposite_vector_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_ORTHOGONAL_VECTOR_3)
CGAL_Kernel_cons(Construct_orthogonal_vector_3,
		 construct_orthogonal_vector_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_PERPENDICULAR_DIRECTION_2)
CGAL_Kernel_cons(Construct_perpendicular_direction_2,
		 construct_perpendicular_direction_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_PERPENDICULAR_LINE_2)
CGAL_Kernel_cons(Construct_perpendicular_line_2,
		 construct_perpendicular_line_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_PERPENDICULAR_LINE_3)
CGAL_Kernel_cons(Construct_perpendicular_line_3,
		 construct_perpendicular_line_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_PERPENDICULAR_PLANE_3)
CGAL_Kernel_cons(Construct_perpendicular_plane_3,
		 construct_perpendicular_plane_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_PERPENDICULAR_VECTOR_2)
CGAL_Kernel_cons(Construct_perpendicular_vector_2,
		 construct_perpendicular_vector_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_PLANE_3)
CGAL_Kernel_cons(Construct_plane_3,
		 construct_plane_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_POINT_ON_2)
CGAL_Kernel_cons(Construct_point_on_2,
		 construct_point_on_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_POINT_ON_3)
CGAL_Kernel_cons(Construct_point_on_3,
		 construct_point_on_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_POINT_2)
CGAL_Kernel_cons(Construct_point_2,
		 construct_point_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_POINT_3)
CGAL_Kernel_cons(Construct_point_3,
		 construct_point_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_PROJECTED_POINT_2)
CGAL_Kernel_cons(Construct_projected_point_2,
		 construct_projected_point_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_PROJECTED_POINT_3)
CGAL_Kernel_cons(Construct_projected_point_3,
		 construct_projected_point_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_PROJECTED_XY_POINT_2)
CGAL_Kernel_cons(Construct_projected_xy_point_2,
		 construct_projected_xy_point_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_RAY_2)
CGAL_Kernel_cons(Construct_ray_2,
		 construct_ray_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_RAY_3)
CGAL_Kernel_cons(Construct_ray_3,
		 construct_ray_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_SCALED_VECTOR_2)
CGAL_Kernel_cons(Construct_scaled_vector_2,
		 construct_scaled_vector_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_SCALED_VECTOR_3)
CGAL_Kernel_cons(Construct_scaled_vector_3,
		 construct_scaled_vector_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_SEGMENT_2)
CGAL_Kernel_cons(Construct_segment_2,
		 construct_segment_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_SEGMENT_3)
CGAL_Kernel_cons(Construct_segment_3,
		 construct_segment_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_SPHERE_3)
CGAL_Kernel_cons(Construct_sphere_3,
		 construct_sphere_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_SUPPORTING_LINE_2)
CGAL_Kernel_cons(Construct_supporting_line_2,
		 construct_supporting_line_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_SUPPORTING_LINE_3)
CGAL_Kernel_cons(Construct_supporting_line_3,
		 construct_supporting_line_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_SUPPORTING_PLANE_3)
CGAL_Kernel_cons(Construct_supporting_plane_3,
		 construct_supporting_plane_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_TETRAHEDRON_3)
CGAL_Kernel_cons(Construct_tetrahedron_3,
		 construct_tetrahedron_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_TRANSLATED_POINT_2)
CGAL_Kernel_cons(Construct_translated_point_2,
		 construct_translated_point_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_TRANSLATED_POINT_3)
CGAL_Kernel_cons(Construct_translated_point_3,
		 construct_translated_point_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_TRIANGLE_2)
CGAL_Kernel_cons(Construct_triangle_2,
		 construct_triangle_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_CONSTRUCT_TRIANGLE_3)
CGAL_Kernel_cons(Construct_triangle_3,
		 construct_triangle_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CONSTRUCT_VECTOR_2)
CGAL_Kernel_cons(Construct_vector_2,
		 construct_vector_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CONSTRUCT_VECTOR_3)
CGAL_Kernel_cons(Construct_vector_3,
		 construct_vector_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CONSTRUCT_VERTEX_2)
CGAL_Kernel_cons(Construct_vertex_2,
		 construct_vertex_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CONSTRUCT_VERTEX_3)
CGAL_Kernel_cons(Construct_vertex_3,
		 construct_vertex_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CONSTRUCT_BBOX_2)
CGAL_Kernel_cons(Construct_bbox_2,
		 construct_bbox_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_CONSTRUCT_BBOX_3)
CGAL_Kernel_cons(Construct_bbox_3,
		 construct_bbox_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COPLANAR_ORIENTATION_3)
CGAL_Kernel_pred(Coplanar_orientation_3,
		 coplanar_orientation_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COPLANAR_SIDE_OF_BOUNDED_CIRCLE_3)
CGAL_Kernel_pred(Coplanar_side_of_bounded_circle_3,
		 coplanar_side_of_bounded_circle_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_COPLANAR_3)
CGAL_Kernel_pred(Coplanar_3,
		 coplanar_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_COUNTERCLOCKWISE_IN_BETWEEN_2)
CGAL_Kernel_pred(Counterclockwise_in_between_2,
		 counterclockwise_in_between_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_DO_INTERSECT_2)
CGAL_Kernel_pred(Do_intersect_2,
		 do_intersect_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_DO_INTERSECT_3)
CGAL_Kernel_cons(Do_intersect_3,
		 do_intersect_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_EQUAL_XY_3)
CGAL_Kernel_pred(Equal_xy_3,
		 equal_xy_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_EQUAL_X_2)
CGAL_Kernel_pred(Equal_x_2,
		 equal_x_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_EQUAL_X_3)
CGAL_Kernel_pred(Equal_x_3,
		 equal_x_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_EQUAL_Y_2)
CGAL_Kernel_pred(Equal_y_2,
		 equal_y_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_EQUAL_Y_3)
CGAL_Kernel_pred(Equal_y_3,
		 equal_y_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_EQUAL_Z_3)
CGAL_Kernel_pred(Equal_z_3,
		 equal_z_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_EQUAL_2)
CGAL_Kernel_pred(Equal_2,
		 equal_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_EQUAL_3)
CGAL_Kernel_pred(Equal_3,
		 equal_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_HAS_ON_BOUNDARY_2)
CGAL_Kernel_pred(Has_on_boundary_2,
		 has_on_boundary_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_HAS_ON_BOUNDARY_3)
CGAL_Kernel_pred(Has_on_boundary_3,
		 has_on_boundary_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_HAS_ON_BOUNDED_SIDE_2)
CGAL_Kernel_pred(Has_on_bounded_side_2,
		 has_on_bounded_side_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_HAS_ON_BOUNDED_SIDE_3)
CGAL_Kernel_pred(Has_on_bounded_side_3,
		 has_on_bounded_side_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_HAS_ON_NEGATIVE_SIDE_2)
CGAL_Kernel_pred(Has_on_negative_side_2,
		 has_on_negative_side_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_HAS_ON_NEGATIVE_SIDE_3)
CGAL_Kernel_pred(Has_on_negative_side_3,
		 has_on_negative_side_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_HAS_ON_POSITIVE_SIDE_2)
CGAL_Kernel_pred(Has_on_positive_side_2,
		 has_on_positive_side_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_HAS_ON_POSITIVE_SIDE_3)
CGAL_Kernel_pred(Has_on_positive_side_3,
		 has_on_positive_side_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_HAS_ON_UNBOUNDED_SIDE_2)
CGAL_Kernel_pred(Has_on_unbounded_side_2,
		 has_on_unbounded_side_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_HAS_ON_UNBOUNDED_SIDE_3)
CGAL_Kernel_pred(Has_on_unbounded_side_3,
		 has_on_unbounded_side_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_HAS_ON_2)
CGAL_Kernel_pred(Has_on_2,
		 has_on_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_HAS_ON_3)
CGAL_Kernel_pred(Has_on_3,
		 has_on_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_INTERSECT_2)
CGAL_Kernel_cons(Intersect_2,
		 intersect_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_INTERSECT_3)
CGAL_Kernel_cons(Intersect_3,
		 intersect_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_IS_DEGENERATE_2)
CGAL_Kernel_pred(Is_degenerate_2,
		 is_degenerate_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_IS_DEGENERATE_3)
CGAL_Kernel_pred(Is_degenerate_3,
		 is_degenerate_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_IS_HORIZONTAL_2)
CGAL_Kernel_pred(Is_horizontal_2,
		 is_horizontal_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_IS_VERTICAL_2)
CGAL_Kernel_pred(Is_vertical_2,
		 is_vertical_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LEFT_TURN_2)
CGAL_Kernel_pred(Left_turn_2,
                 left_turn_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_LESS_DISTANCE_TO_POINT_2)
CGAL_Kernel_pred(Less_distance_to_point_2,
                 less_distance_to_point_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_LESS_DISTANCE_TO_POINT_3)
CGAL_Kernel_pred(Less_distance_to_point_3,
                 less_distance_to_point_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_LESS_ROTATE_CCW_2)
CGAL_Kernel_pred(Less_rotate_ccw_2,
		 less_rotate_ccw_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_LESS_SIGNED_DISTANCE_TO_LINE_2)
CGAL_Kernel_pred(Less_signed_distance_to_line_2,
                 less_signed_distance_to_line_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_LESS_SIGNED_DISTANCE_TO_PLANE_3)
CGAL_Kernel_pred(Less_signed_distance_to_plane_3,
		 less_signed_distance_to_plane_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LESS_XYZ_3)
CGAL_Kernel_pred(Less_xyz_3,
		 less_xyz_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LESS_XY_2)
CGAL_Kernel_pred(Less_xy_2,
		 less_xy_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LESS_XY_3)
CGAL_Kernel_pred(Less_xy_3,
		 less_xy_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LESS_X_2)
CGAL_Kernel_pred(Less_x_2,
		 less_x_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LESS_X_3)
CGAL_Kernel_pred(Less_x_3,
		 less_x_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LESS_YX_2)
CGAL_Kernel_pred(Less_yx_2,
		 less_yx_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LESS_Y_2)
CGAL_Kernel_pred(Less_y_2,
		 less_y_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LESS_Y_3)
CGAL_Kernel_pred(Less_y_3,
		 less_y_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_LESS_Z_3)
CGAL_Kernel_pred(Less_z_3,
		 less_z_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ORIENTATION_2)
CGAL_Kernel_pred(Orientation_2,
		 orientation_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ORIENTATION_3)
CGAL_Kernel_pred(Orientation_3,
		 orientation_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ORIENTED_SIDE_2)
CGAL_Kernel_pred(Oriented_side_2,
		 oriented_side_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || defined(CGAL_CA_ORIENTED_SIDE_3)
CGAL_Kernel_pred(Oriented_side_3,
		 oriented_side_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_SIDE_OF_BOUNDED_CIRCLE_2)
CGAL_Kernel_pred(Side_of_bounded_circle_2,
		 side_of_bounded_circle_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_SIDE_OF_BOUNDED_SPHERE_3)
CGAL_Kernel_pred(Side_of_bounded_sphere_3,
		 side_of_bounded_sphere_3_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_SIDE_OF_ORIENTED_CIRCLE_2)
CGAL_Kernel_pred(Side_of_oriented_circle_2,
		 side_of_oriented_circle_2_object)
#endif

#if !defined(CGAL_CA_LIMITED_INTERFACE) || \
    defined(CGAL_CA_SIDE_OF_ORIENTED_SPHERE_3)
CGAL_Kernel_pred(Side_of_oriented_sphere_3,
		 side_of_oriented_sphere_3_object)
#endif

// deprecated functors were removed ...

#undef CGAL_Kernel_pred
#undef CGAL_Kernel_cons
