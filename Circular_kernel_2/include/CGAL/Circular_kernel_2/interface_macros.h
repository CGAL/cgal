// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)



// This file is intentionally not protected against re-inclusion.
// It's aimed at being included from within a kernel traits class, this
// way we share more code.

// It is the responsability of the including file to correctly set the 2
// macros CGAL_Circular_Kernel_pred and CGAL_Circular_Kernel_cons.
// And they are #undefed at the end of this file.

  CGAL_Circular_Kernel_cons(Compute_squared_radius_2,
  compute_squared_radius_2_object)
  CGAL_Circular_Kernel_cons(Construct_center_2,
  construct_center_2_object)
  CGAL_Circular_Kernel_cons(Get_equation,
  get_equation_object)
  CGAL_Circular_Kernel_cons(Construct_circle_2,
  construct_circle_2_object)
  CGAL_Circular_Kernel_pred(Compare_x_2,
  compare_x_2_object)  
  CGAL_Circular_Kernel_pred(Compare_y_2,
  compare_y_2_object)    
  CGAL_Circular_Kernel_pred(Compare_xy_2,
  compare_xy_2_object)  
  CGAL_Circular_Kernel_pred(Compare_y_at_x_2,
  compare_y_at_x_2_object)  
  CGAL_Circular_Kernel_pred(Compare_y_to_right_2,
  compare_y_to_right_2_object)  
  CGAL_Circular_Kernel_pred(Do_overlap_2,
  do_overlap_2_object)  
  CGAL_Circular_Kernel_pred(Equal_2,
  equal_2_object)
  CGAL_Circular_Kernel_pred(In_x_range_2,
  in_x_range_2_object)
  CGAL_Circular_Kernel_cons(Make_x_monotone_2,
  make_x_monotone_2_object)
  CGAL_Circular_Kernel_cons(Make_xy_monotone_2,
  make_xy_monotone_2_object)
  CGAL_Circular_Kernel_cons(Intersect_2,
  intersect_2_object)
  CGAL_Circular_Kernel_cons(Split_2,
  split_2_object)
  CGAL_Circular_Kernel_cons(Construct_circular_arc_2,
  construct_circular_arc_2_object)
  CGAL_Circular_Kernel_cons(Construct_line_arc_2,
  construct_line_arc_2_object)
  CGAL_Circular_Kernel_cons(Construct_circular_arc_point_2,
  construct_circular_arc_point_2_object)
  CGAL_Circular_Kernel_cons(Compute_circular_x_2,
  compute_circular_x_2_object)  
  CGAL_Circular_Kernel_cons(Compute_circular_y_2,
  compute_circular_y_2_object)
  CGAL_Circular_Kernel_cons(Construct_circular_min_vertex_2,
  construct_circular_min_vertex_2_object)
  CGAL_Circular_Kernel_cons(Construct_circular_max_vertex_2,
  construct_circular_max_vertex_2_object)
  CGAL_Circular_Kernel_cons(Construct_circular_source_vertex_2,
  construct_circular_source_vertex_2_object)
  CGAL_Circular_Kernel_cons(Construct_circular_target_vertex_2,
  construct_circular_target_vertex_2_object)
  CGAL_Circular_Kernel_pred(Is_x_monotone_2,
  is_x_monotone_2_object)  
  CGAL_Circular_Kernel_pred(Is_y_monotone_2,
  is_y_monotone_2_object)
  CGAL_Circular_Kernel_pred(Is_vertical_2,
  is_vertical_2_object) 
  CGAL_Circular_Kernel_pred(Has_on_2,
  has_on_2_object)
  
  CGAL_Circular_Kernel_pred(Has_on_bounded_side_2, has_on_bounded_side_2_object)
  CGAL_Circular_Kernel_pred(Has_on_unbounded_side_2, has_on_unbounded_side_2_object)
  CGAL_Circular_Kernel_pred(Bounded_side_2, bounded_side_2_object)
  CGAL_Circular_Kernel_pred(Do_intersect_2, do_intersect_2_object)

#ifndef CGAL_NO_DEPRECATED_CODE
	CGAL_Circular_Kernel_cons(Construct_supporting_circle_2,
	  construct_supporting_circle_2_object)
	CGAL_Circular_Kernel_cons(Construct_supporting_line_2,
	  construct_supporting_line_2_object)
#endif
	
  CGAL_Circular_Kernel_cons(Construct_bbox_2,
  construct_bbox_2_object)


#undef CGAL_Circular_Kernel_pred
#undef CGAL_Circular_Kernel_cons
