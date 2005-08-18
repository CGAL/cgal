// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Circular_arc_traits.h

#ifndef CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_H
#define CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_H

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/global_functions_on_circular_arcs_2.h>

namespace CGAL {

/// Traits class for CGAL::Arrangement_2 (and similar) based on a CurvedKernel.

template < typename CurvedKernel >
class Circular_arc_traits {

  CurvedKernel ck;

public:

  typedef CurvedKernel Kernel;
  typedef typename CurvedKernel::Circular_arc_2  Curve_2;
  typedef typename CurvedKernel::Circular_arc_2  X_monotone_curve_2;

  typedef typename CurvedKernel::Circular_arc_endpoint_2      Point;
  typedef typename CurvedKernel::Circular_arc_endpoint_2      Point_2;

  typedef CGAL::Tag_false                        Has_left_category;
  typedef CGAL::Tag_false 			 Has_merge_category;

  Circular_arc_traits(const CurvedKernel &k = CurvedKernel())
    : ck(k) {}

  typedef typename CurvedKernel::Compare_x_2            Compare_x_2;
  typedef typename CurvedKernel::Compare_xy_2           Compare_xy_2;
  typedef typename CurvedKernel::Compare_y_at_x_2       Compare_y_at_x_2;
  typedef typename CurvedKernel::Compare_y_to_right_2   Compare_y_at_x_right_2;
  typedef typename CurvedKernel::Construct_max_vertex_2 Construct_max_vertex_2;
  typedef typename CurvedKernel::Construct_min_vertex_2 Construct_min_vertex_2;     
  typedef typename CurvedKernel::Equal_2                Equal_2;
  typedef typename CurvedKernel::Make_x_monotone_2      Make_x_monotone_2;
  typedef typename CurvedKernel::Split_2                Split_2;
  typedef typename CurvedKernel::Construct_intersections_2 Intersect_2;


  class Is_vertical_2
  {
  public:
    typedef bool result_type;

    // TO BE IMPLEMENTED !!!!!!!
    bool operator() (const X_monotone_curve_2& cv) const
    {
      return false;
    }
  };

  
  Compare_x_2 compare_x_2_object() const
  { return ck.compare_x_2_object(); }

  Compare_xy_2 compare_xy_2_object() const
  { return ck.compare_xy_2_object(); }

  Compare_y_at_x_2 compare_y_at_x_2_object() const 
  { return ck.compare_y_at_x_2_object(); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const 
  { return ck.compare_y_to_right_2_object(); }

  Equal_2 equal_2_object() const
  { return ck.equal_2_object(); }

  Make_x_monotone_2 make_x_monotone_2_object() const
  { return ck.make_x_monotone_2_object(); }

  Split_2 split_2_object() const
  { return ck.split_2_object(); }

  Intersect_2 intersect_2_object() const
    { return ck.construct_intersections_2_object(); }
    
  Construct_max_vertex_2 construct_max_vertex_2_object() const
    { return ck.construct_max_vertex_2_object(); }
       
  Construct_min_vertex_2 construct_min_vertex_2_object() const
    { return ck.construct_min_vertex_2_object(); }
       
  Is_vertical_2 is_vertical_2_object() const
    { return Is_vertical_2();}


};

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_H
