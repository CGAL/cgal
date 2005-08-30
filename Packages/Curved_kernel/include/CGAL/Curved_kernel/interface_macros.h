// Copyright (c) 2000-2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Herve Bronnimann, Sylvain Pion, Susan Hert

// This file is intentionally not protected against re-inclusion.
// It's aimed at being included from within a kernel traits class, this
// way we share more code.

// It is the responsability of the including file to correctly set the 2
// macros CGAL_Kernel_pred and CGAL_Kernel_cons.
// And they are #undefed at the end of this file.

  CGAL_Curved_Kernel_cons(Get_equation,
  get_equation_object)
  CGAL_Curved_Kernel_cons(Construct_circle_2,
  construct_circle_2_object)
  CGAL_Curved_Kernel_pred(Compare_x_2,
  compare_x_2_object)  
  CGAL_Curved_Kernel_pred(Compare_y_2,
  compare_y_2_object)    
  CGAL_Curved_Kernel_pred(Compare_xy_2,
  compare_xy_2_object)  
  CGAL_Curved_Kernel_pred(Compare_y_at_x_2,
  compare_y_at_x_2_object)  
  CGAL_Curved_Kernel_pred(Compare_y_to_right_2,
  compare_y_to_right_2_object)  
  CGAL_Curved_Kernel_pred(Do_overlap_2,
  do_overlap_2_object)  
  CGAL_Curved_Kernel_pred(Equal_2,
  equal_2_object)
  CGAL_Curved_Kernel_pred(In_range_2,
  in_range_2_object)
  CGAL_Curved_Kernel_cons(Make_x_monotone_2,
  make_x_monotone_2_object)
  CGAL_Curved_Kernel_cons(Intersect_2,
  intersect_2_object)
  CGAL_Curved_Kernel_cons(Split_2,
  split_2_object)
  CGAL_Curved_Kernel_cons(Construct_circular_arc_2,
  construct_circular_arc_2_object)
  CGAL_Curved_Kernel_cons(Construct_line_arc_2,
  construct_line_arc_2_object)
  CGAL_Curved_Kernel_cons(Construct_circular_arc_point_2,
  construct_circular_arc_point_2_object)
  CGAL_Curved_Kernel_cons(Compute_x_2,
  compute_x_2_object)  
  CGAL_Curved_Kernel_cons(Compute_y_2,
  compute_y_2_object)
  CGAL_Curved_Kernel_cons(Construct_min_vertex_2,
  construct_min_vertex_2_object)
  CGAL_Curved_Kernel_cons(Construct_max_vertex_2,
  construct_max_vertex_2_object)
  CGAL_Curved_Kernel_cons(Construct_source_vertex_2,
  construct_source_vertex_2_object)
  CGAL_Curved_Kernel_cons(Construct_target_vertex_2,
  construct_target_vertex_2_object)
  CGAL_Curved_Kernel_pred(Is_x_monotone_2,
  is_x_monotone_2_object)  
  CGAL_Curved_Kernel_pred(Is_y_monotone_2,
  is_y_monotone_2_object)
  CGAL_Curved_Kernel_pred(Is_vertical_2,
  is_vertical_2_object) 
  CGAL_Curved_Kernel_cons(Construct_supporting_circle_2,
  construct_supporting_circle_2_object)
  CGAL_Curved_Kernel_cons(Construct_supporting_line_2,
  construct_supporting_line_2_object)
  CGAL_Curved_Kernel_cons(Construct_bbox_2,
  construct_bbox_2_object)


#undef CGAL_Curved_Kernel_pred
#undef CGAL_Curved_Kernel_cons
