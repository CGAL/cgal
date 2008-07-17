// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//             Pedro Machado    <tashimir@gmail.com>
//             Julien Hazebrouck
//             Damien Leroy

// This file is intentionally not protected against re-inclusion.
// It's aimed at being included from within a kernel traits class, this
// way we share more code.

// It is the responsability of the including file to correctly set the 2
// macros CGAL_Kernel_pred and CGAL_Kernel_cons.
// And they are #undefed at the end of this file.

  
  //TAG_SEB
  CGAL_Spherical_Kernel_cons(Compute_type_of_circle_on_reference_sphere_3,Compute_type_of_circle_on_reference_sphere_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circle_center_coefficient_3,compute_circle_center_coefficient_3_object)
  CGAL_Spherical_Kernel_cons(Compute_extremal_point_z,compute_extremal_point_z_object)
  CGAL_Spherical_Kernel_cons(Compute_supporting_sphere_radius_3,compute_supporting_sphere_radius_3_object)
  CGAL_Spherical_Kernel_cons(Compute_supporting_sphere_squared_radius_3,compute_supporting_sphere_squared_radius_3_object)
  CGAL_Spherical_Kernel_cons(Compute_supporting_sphere_center_3,compute_supporting_sphere_center_3_object)
  CGAL_Spherical_Kernel_cons(Compute_reference_sphere_3,compute_reference_sphere_3_object)  
  CGAL_Spherical_Kernel_cons(Compute_circular_theta_rep_3, compute_circular_theta_rep_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_theta_3, compute_circular_theta_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_hq_3, compute_circular_hq_3_object)
  CGAL_Spherical_Kernel_cons(Compute_theta_hq_3, compute_theta_hq_3_object)
  CGAL_Spherical_Kernel_cons(Compute_theta_ftheta_3, compute_theta_ftheta_3_object)
  CGAL_Spherical_Kernel_cons(Compute_radius_sphere_with_radius_3, construct_radius_sphere_with_radius_3_object)
  CGAL_Spherical_Kernel_cons(Compute_supporting_circle_on_reference_sphere_3, compute_supporting_circle_on_reference_sphere_3_object)
  CGAL_Spherical_Kernel_cons(Compute_half_circle_position_3,compute_half_circle_position_3)
  CGAL_Spherical_Kernel_cons(Construct_half_circle_on_reference_sphere_3, construct_half_circle_on_reference_sphere_3_object)
  CGAL_Spherical_Kernel_cons(Construct_sphere_with_radius_3, construct_sphere_with_radius_3_object)
  CGAL_Spherical_Kernel_cons(Construct_theta_rep, construct_theta_rep_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_arc_point_on_reference_sphere_3, construct_circular_arc_point_on_reference_sphere_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circle_on_reference_sphere_3,construct_circle_on_reference_sphere_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_arc_on_reference_sphere_3,construct_circular_arc_on_reference_sphere_3_object)
  CGAL_Spherical_Kernel_cons(Theta_extremal_point_3,theta_extremal_point_3_object)
  CGAL_Spherical_Kernel_cons(Compare_theta_3,compare_theta_3_object)
  CGAL_Spherical_Kernel_cons(Compare_theta_z_3,compare_theta_z_3_object)
  CGAL_Spherical_Kernel_cons(Compare_z_at_theta_3,compare_z_at_theta_3_object)
  CGAL_Spherical_Kernel_cons(Compare_z_to_left_3,compare_z_to_left_3_object)

  CGAL_Spherical_Kernel_cons(Get_equation, get_equation_object) 
  CGAL_Spherical_Kernel_cons(Construct_circular_arc_point_3, construct_circular_arc_point_3_object)
  CGAL_Spherical_Kernel_cons(Construct_sphere_3, construct_sphere_3_object)
  CGAL_Spherical_Kernel_cons(Construct_plane_3, construct_plane_3_object)
  CGAL_Spherical_Kernel_cons(Construct_line_3, construct_line_3_object)
  CGAL_Spherical_Kernel_cons(Construct_line_arc_3, construct_line_arc_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_arc_3, construct_circular_arc_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circle_3, construct_circle_3_object)
  CGAL_Spherical_Kernel_cons(Construct_diametral_sphere_3, construct_diametral_sphere_3_object)
  CGAL_Spherical_Kernel_cons(Construct_supporting_plane_3, construct_supporting_plane_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_x_3, compute_circular_x_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_y_3, compute_circular_y_3_object)
  CGAL_Spherical_Kernel_cons(Compute_circular_z_3, compute_circular_z_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_min_vertex_3, construct_circular_min_vertex_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_max_vertex_3, construct_circular_max_vertex_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_source_vertex_3, construct_circular_source_vertex_3_object)
  CGAL_Spherical_Kernel_cons(Construct_circular_target_vertex_3, construct_circular_target_vertex_3_object)
  CGAL_Spherical_Kernel_cons(Construct_supporting_line_3, construct_supporting_line_3_object)
  CGAL_Spherical_Kernel_cons(Construct_supporting_circle_3, construct_supporting_circle_3_object)
  CGAL_Spherical_Kernel_cons(Intersect_3, intersect_3_object)
  CGAL_Spherical_Kernel_cons(Construct_bbox_3, construct_bbox_3_object)
  CGAL_Spherical_Kernel_cons(Split_3, split_3_object)
  CGAL_Spherical_Kernel_cons(Compute_approximate_squared_length_3, compute_approximate_squared_length_3_object)
  CGAL_Spherical_Kernel_cons(Compute_approximate_angle_3, compute_approximate_angle_3_object)

  CGAL_Spherical_Kernel_pred(Compare_x_3, compare_x_3_object)
  CGAL_Spherical_Kernel_pred(Compare_y_3, compare_y_3_object)
  CGAL_Spherical_Kernel_pred(Compare_z_3, compare_z_3_object)
  CGAL_Spherical_Kernel_pred(Compare_xy_3, compare_xy_3_object)
  CGAL_Spherical_Kernel_pred(Compare_xyz_3, compare_xyz_3_object)
  CGAL_Spherical_Kernel_pred(Equal_3, equal_3_object)
  CGAL_Spherical_Kernel_pred(Has_on_3, has_on_3_object)
  CGAL_Spherical_Kernel_pred(Has_on_bounded_side_3, has_on_bounded_side_3_object)
  CGAL_Spherical_Kernel_pred(Has_on_unbounded_side_3, has_on_unbounded_side_3_object)
  CGAL_Spherical_Kernel_pred(Bounded_side_3, bounded_side_3_object)
  CGAL_Spherical_Kernel_pred(Do_overlap_3, do_overlap_3_object)

#undef CGAL_Spherical_Kernel_pred
#undef CGAL_Spherical_Kernel_cons
