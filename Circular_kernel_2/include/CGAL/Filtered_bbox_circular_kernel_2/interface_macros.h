// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://pmachado@scm.gforge.inria.fr/svn/cgal/trunk/Circular_kernel_2/include/CGAL/Circular_kernel_2/interface_macros.h $
// $Id: interface_macros.h 44444 2008-07-25 15:02:26Z pmachado $
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
// macros CGAL_Filtered_Bbox_Circular_Kernel_pred and CGAL_Filtered_Bbox_Circular_Kernel_cons.
// And they are #undefed at the end of this file.

  CGAL_Filtered_Bbox_Circular_Kernel_pred(Compare_x_2, compare_x_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_pred(Compare_y_2, compare_y_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_pred(Compare_xy_2, compare_xy_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_pred(Has_on_2, has_on_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_pred(Compare_y_at_x_2, compare_y_at_x_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_pred(Compare_y_to_right_2, compare_y_to_right_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_pred(Do_overlap_2, do_overlap_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_pred(Equal_2, equal_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_pred(In_x_range_2, in_x_range_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_pred(Do_intersect_2, do_intersect_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_pred(Is_vertical_2, is_vertical_2_object)

  CGAL_Filtered_Bbox_Circular_Kernel_cons(Construct_line_arc_2,
  construct_line_arc_2_object_object)
  CGAL_Filtered_Bbox_Circular_Kernel_cons(Construct_circular_arc_2,
  construct_circular_arc_2_object_object)
  CGAL_Filtered_Bbox_Circular_Kernel_cons(Construct_circular_source_vertex_2,
  construct_circular_source_vertex_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_cons(Construct_circular_target_vertex_2,
  construct_circular_target_vertex_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_cons(Construct_circular_min_vertex_2,
  construct_circular_min_vertex_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_cons(Construct_circular_max_vertex_2,
  construct_circular_max_vertex_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_cons(Make_x_monotone_2,
  make_x_monotone_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_cons(Intersect_2,
  intersect_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_cons(Split_2,
  split_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_cons(Construct_bbox_2,
  construct_bbox_2_object)
  CGAL_Filtered_Bbox_Circular_Kernel_cons(Construct_circular_arc_point_2,
  construct_circular_arc_point_2_object)

#undef CGAL_Filtered_Bbox_Circular_Kernel_pred
#undef CGAL_Filtered_Bbox_Circular_Kernel_cons

