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

// file : include/CGAL/Circular_kernel.h

#ifndef CGAL_CIRCULAR_KERNEL_H
#define CGAL_CIRCULAR_KERNEL_H

#include <CGAL/Curved_kernel/Circular_arc_endpoint_2.h>
#include <CGAL/Curved_kernel/Circular_arc_2.h>
#include <CGAL/Circular_arc_2.h>
#include <CGAL/Circular_arc_endpoint_2.h>

#include <CGAL/Curved_kernel/function_objects_on_circle_2.h>
#include <CGAL/global_functions_on_circle_2.h>

#include <CGAL/Curved_kernel/function_objects_polynomial_circular.h>
#include <CGAL/global_functions_on_circular_arcs_2.h>

namespace CGAL {
namespace CGALi {

template < class CurvedKernel >
struct Curved_kernel_base
// takes classes in internal sub-namespace
{
  typedef CGALi::Circular_arc_2<CurvedKernel>               Circular_arc_2;
  typedef CGALi::Circular_arc_endpoint_2<CurvedKernel>      Circular_arc_endpoint_2;
};

} // namespace CGALi

template < class LinearKernel, class AlgebraicKernel >
struct Curved_kernel
  : public LinearKernel,
    // there should be a derivation from
    // LinearKernel::Kernel_base<Self> to have types equalities for
    // the Linearkernel types
    public CGALi::Curved_kernel_base<Curved_kernel<LinearKernel,AlgebraicKernel> >
{
  typedef LinearKernel                                    Linear_kernel;
  typedef AlgebraicKernel                                 Algebraic_kernel;

  typedef Curved_kernel<LinearKernel,AlgebraicKernel>      Self;
  typedef CGALi::Curved_kernel_base<Self>                  Kernel_base;

  typedef typename LinearKernel::RT                       RT;
  typedef typename LinearKernel::FT                       FT;

  typedef typename Algebraic_kernel::Root_of_2            Root_of_2;
  typedef typename Algebraic_kernel::Polynomial_for_circles_2_2
                                                  Polynomial_for_circles_2_2;
  typedef typename Algebraic_kernel::Polynomial_1_2
                                                  Polynomial_1_2;

  // public classes

  typedef typename Linear_kernel::Line_2                  Line_2;
  typedef typename Linear_kernel::Circle_2                Circle_2;
  typedef typename Linear_kernel::Conic_2                 Conic_2;
  typedef typename Linear_kernel::Point_2                 Point_2;

  typedef CGAL::Circular_arc_2<Self>               Circular_arc_2;
  typedef CGAL::Circular_arc_endpoint_2<Self>      Circular_arc_endpoint_2;

  typedef CircularFunctors::Construct_circle_2<Self>   Construct_circle_2;
  typedef CircularFunctors::Get_equation<Self>         Get_equation;

  // Function objects
  typedef CircularFunctors::Compare_x_2<Self>
                              Compare_x_2;
  typedef CircularFunctors::Compare_y_2<Self>
                              Compare_y_2;
  typedef CircularFunctors::Compare_xy_2<Self>
                              Compare_xy_2;
  typedef CircularFunctors::Compare_y_at_x_2<Self>
                              Compare_y_at_x_2;
  typedef CircularFunctors::Compare_y_to_right_2<Self>
                              Compare_y_to_right_2;
  typedef CircularFunctors::Do_overlap_2<Self>
                              Do_overlap_2;
  typedef CircularFunctors::Equal_2<Self>
                              Equal_2;
  typedef CircularFunctors::In_range_2<Self>
                              In_range_2;
  typedef CircularFunctors::Make_x_monotone_2<Self>
                              Make_x_monotone_2;
  typedef CircularFunctors::Construct_intersections_2<Self>
                              Construct_intersections_2;
  typedef CircularFunctors::Split_2<Self>
                              Split_2;
  typedef CircularFunctors::Nearest_intersection_to_right_2<Self>
                              Nearest_intersection_to_right_2;

  Get_equation
  get_equation_object() const
    { return Get_equation(); }

  Construct_circle_2
  construct_circle_2_object() const
    { return Construct_circle_2(); }

  Compare_x_2
  compare_x_2_object() const
  { return Compare_x_2(); }

  Compare_y_2
  compare_y_2_object() const
  { return Compare_y_2(); }

  Compare_xy_2
  compare_xy_2_object() const
  { return Compare_xy_2(); }

  Compare_y_at_x_2
  compare_y_at_x_2_object() const 
  { return Compare_y_at_x_2(); }

  Compare_y_to_right_2
  compare_y_to_right_2_object() const
  { return Compare_y_to_right_2(); }

  Do_overlap_2
  do_overlap_2_object() const
  { return Do_overlap_2(); }

  Equal_2
  equal_2_object() const
  { return Equal_2(); }

  In_range_2
  in_range_2_object() const
  { return In_range_2(); }

  Make_x_monotone_2
  make_x_monotone_2_object() const
  { return Make_x_monotone_2(); }

  Construct_intersections_2
  construct_intersections_2_object() const
    { return Construct_intersections_2(); }

  Nearest_intersection_to_right_2
  nearest_intersection_to_right_2_object() const
  { return Nearest_intersection_to_right_2(); }

  Split_2
  split_2_object() const
  { return Split_2(); }

};

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_H
