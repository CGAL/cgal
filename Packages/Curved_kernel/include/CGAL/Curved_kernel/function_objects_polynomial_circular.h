// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
//                     National & Kapodistrian University of Athens (Greece).
// All rights reserved.
//
// Authors : Monique Teillaud     <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion         <Sylvain.Pion@sophia.inria.fr>
//           Athanasios Kakargias <grad0460@di.uoa.gr>
// 
// Partially supported by INRIA's project "CALAMATA", a
// bilateral collaboration with National Kapodistrian University of
// Athens.
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file: include/CGAL/Curved_kernel/function_objects_polynomial.h

#ifndef CGAL_CURVED_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_CIRCULAR_H
#define CGAL_CURVED_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_CIRCULAR_H

#include <CGAL/Curved_kernel/internal_functions_on_circular_arc_2.h>

namespace CGAL {
namespace CircularFunctors {

  template < class CK >
  class Compare_x_2
    : public CK::Linear_kernel::Compare_x_2
  {
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return compare_x<CK>(p0, p1); }
  };

  template < class CK >
  class Compare_y_2
    : public CK::Linear_kernel::Compare_y_2
  {
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return compare_y<CK>(p0, p1); }

  };

  template < class CK >
  class Compare_xy_2
    : public CK::Linear_kernel::Compare_xy_2
  {
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return compare_xy<CK>(p0, p1); }

  };

  template < class CK >
  class In_range_2
  {
    typedef typename CK::Circular_arc_2          Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;

  public:
    typedef bool result_type;

    result_type
    operator()(const Circular_arc_2 &a, const Circular_arc_endpoint_2 &p) const
    { return point_in_range<CK>(a, p); }

  };

  template < class CK >
  class Compare_y_to_right_2
  {
    typedef typename CK::Circular_arc_2           Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2  Circular_arc_endpoint_2;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator()(const Circular_arc_2 &a1,
               const Circular_arc_2 &a2,
               const Circular_arc_endpoint_2 &p) const
    { return compare_y_to_right<CK>(a1, a2, p); }

  };

  template < class CK >
  class Equal_2
    : public CK::Linear_kernel::Equal_2
  {
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
    typedef typename CK::Circular_arc_2  Circular_arc_2;

  public:
    typedef bool result_type;

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return equal<CK>(p0, p1); }
    
    result_type
    operator() (const Circular_arc_2 &a0, const Circular_arc_2 &a1) const
    { return equal<CK>(a0, a1); }
    
   };

  template < class CK >
  class Compare_y_at_x_2
  {
    typedef typename CK::Circular_arc_2          Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator() (const Circular_arc_endpoint_2 &p,
                const Circular_arc_2 &A1) const
    { return compare_y_at_x<CK>(p, A1); }

  };

  template < class CK >
  class Do_overlap_2
  {
    typedef typename CK::Circular_arc_2 Circular_arc_2;

  public:
    typedef bool result_type;

    result_type
    operator() (const Circular_arc_2 &A1, const Circular_arc_2 &A2) const
    { return do_overlap<CK>(A1, A2); }
    
  };

  template < class CK >
  class Make_x_monotone_2
  {
    typedef typename CK::Circular_arc_2 Circular_arc_2;

  public:
    // typedef OutputIterator result_type;

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_2 &A, OutputIterator res)
      { return make_x_monotone<CK> (A, res); }

  };

  template < class CK >
  class Construct_intersections_2
  {
    public:

    typedef typename CK::Circle_2                 Circle;
    typedef typename CK::Circular_arc_2           Circular_arc;

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle & c1, const Circle & c2, OutputIterator res)
      { return construct_intersections_2<CK> (c1,c2,res); }

     template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc & c1, const Circular_arc & c2, 
	       OutputIterator res)
      { return construct_intersections_2<CK> (c1,c2,res); }
    
  };

  template < class CK >
  class Nearest_intersection_to_right_2
  {
    typedef typename CK::Circular_arc_2          Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;

  public:
    typedef bool result_type;

    result_type
    operator()(const Circular_arc_2 &A1, const Circular_arc_2 &A2,
               const Circular_arc_endpoint_2 &pt,
               Circular_arc_endpoint_2 &p1, Circular_arc_endpoint_2 &p2) const
    { return nearest_intersection_to_right<CK>(A1, A2, pt, p1, p2); }

  };

  template < class CK >
  class Split_2
  {
    typedef typename CK::Circular_arc_2          Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;

  public:
    typedef void result_type;

    result_type
    operator()(const Circular_arc_2 &A, 
	       const Circular_arc_endpoint_2 &p,
	       Circular_arc_2 &ca1, Circular_arc_2 &ca2) const
    { return split<CK>(A, p, ca1, ca2); }
  };

} // namespace CircularFunctors
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_CIRCULAR_H
