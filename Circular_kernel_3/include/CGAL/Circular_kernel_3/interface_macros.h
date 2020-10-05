// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado,
//             Sebastien Loriot, Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a
// STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

// This file is intentionally not protected against re-inclusion.
// It's aimed at being included from within a kernel traits class, this
// way we share more code.

// It is the responsability of the including file to correctly set the 2
// macros CGAL_Kernel_pred and CGAL_Kernel_cons.
// And they are #undefed at the end of this file.

  CGAL_Spherical_Kernel_cons(Get_equation, get_equation_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_arc_point_3, construct_circular_arc_point_3_object)
  CGAL_Spherical_Kernel_overload_cons(Construct_sphere_3, construct_sphere_3_object)
  CGAL_Spherical_Kernel_overload_cons(Construct_plane_3, construct_plane_3_object)
  CGAL_Spherical_Kernel_cons(Construct_line_arc_3, construct_line_arc_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_arc_3, construct_circular_arc_3_object)
  CGAL_Spherical_Kernel_overload_cons(Construct_circle_3, construct_circle_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_x_3, compute_circular_x_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_y_3, compute_circular_y_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_z_3, compute_circular_z_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_min_vertex_3, construct_circular_min_vertex_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_max_vertex_3, construct_circular_max_vertex_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_source_vertex_3, construct_circular_source_vertex_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_target_vertex_3, construct_circular_target_vertex_3_object)
  CGAL_Spherical_Kernel_overload_cons(Construct_line_3, construct_line_3_object)
  CGAL_Spherical_Kernel_overload_cons(Intersect_3, intersect_3_object)
  CGAL_Spherical_Kernel_overload_cons(Construct_bbox_3, construct_bbox_3_object)
  CGAL_Spherical_Kernel_cons(Split_3, split_3_object)
  CGAL_Spherical_Kernel_overload_cons(Compute_approximate_squared_length_3, compute_approximate_squared_length_3_object)
  CGAL_Spherical_Kernel_cons(Compute_approximate_angle_3, compute_approximate_angle_3_object)

  CGAL_Spherical_Kernel_overload_pred(Compare_x_3, compare_x_3_object)
  CGAL_Spherical_Kernel_overload_pred(Compare_y_3, compare_y_3_object)
  CGAL_Spherical_Kernel_overload_pred(Compare_z_3, compare_z_3_object)
  CGAL_Spherical_Kernel_overload_pred(Compare_xy_3, compare_xy_3_object)
  CGAL_Spherical_Kernel_overload_pred(Compare_xyz_3, compare_xyz_3_object)
  CGAL_Spherical_Kernel_overload_pred(Equal_3, equal_3_object)
  CGAL_Spherical_Kernel_overload_pred(Has_on_3, has_on_3_object)
  CGAL_Spherical_Kernel_overload_pred(Has_on_bounded_side_3, has_on_bounded_side_3_object)
  CGAL_Spherical_Kernel_overload_pred(Has_on_unbounded_side_3, has_on_unbounded_side_3_object)
  CGAL_Spherical_Kernel_overload_pred(Bounded_side_3, bounded_side_3_object)
  CGAL_Spherical_Kernel_pred(Do_overlap_3, do_overlap_3_object)
  CGAL_Spherical_Kernel_overload_pred(Do_intersect_3, do_intersect_3_object)
  CGAL_Spherical_Kernel_pred_on_sphere(Is_theta_monotone_3, is_theta_monotone_3_object)
  CGAL_Spherical_Kernel_pred_on_sphere(Compare_theta_3, compare_theta_3_object)
  CGAL_Spherical_Kernel_pred_on_sphere(Compare_theta_z_3, compare_theta_z_3_object)
  CGAL_Spherical_Kernel_pred_on_sphere(Make_theta_monotone_3, make_theta_monotone_3_object)
  CGAL_Spherical_Kernel_pred_on_sphere(Compare_z_to_right_3, compare_z_to_right_3_object)
  CGAL_Spherical_Kernel_pred_on_sphere(Compare_z_at_theta_3, compare_z_at_theta_3_object)

#undef CGAL_Spherical_Kernel_pred
#undef CGAL_Spherical_Kernel_overload_pred
#undef CGAL_Spherical_Kernel_cons
#undef CGAL_Spherical_Kernel_overload_cons
#undef CGAL_Spherical_Kernel_pred_on_sphere
