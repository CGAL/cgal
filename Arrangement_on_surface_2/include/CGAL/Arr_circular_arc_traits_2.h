// Copyright (c) 2003,2004,2005,2006,2007,2008,2009,2010,2011 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_CIRCULAR_ARC_TRAITS_2_H
#define CGAL_CIRCULAR_KERNEL_CIRCULAR_ARC_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * This file was developed at Inria, France, and copied over to the
 * Arrangement_2 package, which it is now part of. It contains a traits
 * class for the arrangement package that handles circular curves.
 * It is based on the circular kernel.
 *
 * \todo Fix the circular-kernel make-x-monotone functor to use modern variant
 *       instead of the legacy CGAL::Object. Then, eliminate the special
 *       implementation here and directly use the kernel functor instead.
 */

#include <CGAL/basic.h>
#include <CGAL/global_functions_circular_kernel_2.h>
#include <CGAL/Arr_tags.h>

namespace CGAL {

namespace internal{
template <class CircularKernel>
class Non_x_monotonic_Circular_arc_2
  : public CircularKernel::Circular_arc_2
{
  typedef typename CircularKernel::FT             FT;
  typedef typename CircularKernel::Point_2        Point_2;
  typedef typename CircularKernel::Line_2         Line_2;
  typedef typename CircularKernel::Circle_2       Circle_2;
  typedef typename CircularKernel::Circular_arc_point_2
                                                Circular_arc_point_2;

  typedef typename CircularKernel::Circular_arc_2 Base;

public:
  Non_x_monotonic_Circular_arc_2(): Base(){}

  Non_x_monotonic_Circular_arc_2(const Circle_2 &c): Base(c){}
  // Not Documented
  Non_x_monotonic_Circular_arc_2(const Circle_2 &support,
                                 const Line_2 &l1, const bool b_l1,
                                 const Line_2 &l2, const bool b_l2)
    : Base(support,l1,b_l1,l2,b_l2){}

  // Not Documented
  Non_x_monotonic_Circular_arc_2(const Circle_2 &c,
                                 const Circle_2 &c1, const bool b_1,
                                 const Circle_2 &c2, const bool b_2)
    : Base(c,c1,b_1,c2,b_2)
  {}

  Non_x_monotonic_Circular_arc_2(const Point_2 &start,
                                 const Point_2 &middle,
                                 const Point_2 &end)
    : Base(start,middle,end)
  {}

  Non_x_monotonic_Circular_arc_2(const Circle_2 &support,
                                 const Circular_arc_point_2 &begin,
                                 const Circular_arc_point_2 &end)
    : Base(support,begin,end)
  {}

  Non_x_monotonic_Circular_arc_2(const Point_2 &start,
                                 const Point_2 &end,
                                 const FT &bulge)
    : Base(start,end,bulge)
  {}

 Non_x_monotonic_Circular_arc_2(const Base& a) : Base(a) {}
};

} //namespace internal

// Traits class for CGAL::Arrangement_2 (and similar) based on a
// CircularKernel.

template < typename CircularKernel >
class Arr_circular_arc_traits_2 {

  CircularKernel ck;

public:

  typedef CircularKernel Kernel;
  typedef internal::Non_x_monotonic_Circular_arc_2<CircularKernel>  Curve_2;
  typedef typename CircularKernel::Circular_arc_2                   X_monotone_curve_2;

  typedef typename CircularKernel::Circular_arc_point_2 Point;
  typedef typename CircularKernel::Circular_arc_point_2 Point_2;

  typedef unsigned int                           Multiplicity;

  typedef CGAL::Tag_false                        Has_left_category;
  typedef CGAL::Tag_false                        Has_merge_category;
  typedef CGAL::Tag_false                        Has_do_intersect_category;

  typedef Arr_oblivious_side_tag                 Left_side_category;
  typedef Arr_oblivious_side_tag                 Bottom_side_category;
  typedef Arr_oblivious_side_tag                 Top_side_category;
  typedef Arr_oblivious_side_tag                 Right_side_category;

  Arr_circular_arc_traits_2(const CircularKernel &k = CircularKernel())
    : ck(k) {}

  typedef typename CircularKernel::Compare_x_2          Compare_x_2;
  typedef typename CircularKernel::Compare_xy_2         Compare_xy_2;
  typedef typename CircularKernel::Compare_y_at_x_2     Compare_y_at_x_2;
  typedef typename CircularKernel::Compare_y_to_right_2 Compare_y_at_x_right_2;
  typedef typename CircularKernel::Construct_circular_max_vertex_2
                                                        Construct_max_vertex_2;
  typedef typename CircularKernel::Construct_circular_min_vertex_2
                                                        Construct_min_vertex_2;
  typedef typename CircularKernel::Equal_2              Equal_2;
  // typedef typename CircularKernel::Make_x_monotone_2    Make_x_monotone_2;
  typedef typename CircularKernel::Split_2              Split_2;
  typedef typename CircularKernel::Intersect_2          Intersect_2;
  typedef typename CircularKernel::Is_vertical_2        Is_vertical_2;

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

  // Make_x_monotone_2 make_x_monotone_2_object() const
  // { return ck.make_x_monotone_2_object(); }

  Split_2 split_2_object() const
  { return ck.split_2_object(); }

  Intersect_2 intersect_2_object() const
    { return ck.intersect_2_object(); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
    { return ck.construct_circular_max_vertex_2_object(); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
    { return ck.construct_circular_min_vertex_2_object(); }

  Is_vertical_2 is_vertical_2_object() const
    { return ck.is_vertical_2_object();  }


  //! A functor for subdividing curves into x-monotone curves.
  class Make_x_monotone_2 {
  public:
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& arc, OutputIterator oi) const
    {
      typedef boost::variant<Point_2, X_monotone_curve_2>
        Make_x_monotone_result;

      std::vector<CGAL::Object> objs;
      CircularKernel().make_x_monotone_2_object()(arc, std::back_inserter(objs));
      for (const auto& obj : objs) {
        if (const auto* p = CGAL::object_cast<Point_2>(&obj)) {
          *oi++ = Make_x_monotone_result(*p);
          continue;
        }
        if (const auto* xcv = CGAL::object_cast<X_monotone_curve_2>(&obj)) {
          *oi++ = Make_x_monotone_result(*xcv);
          continue;
        }
        CGAL_error();
      }
      return oi;
    }
  };

  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(); }
};

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_CIRCULAR_KERNEL_CIRCULAR_ARC_TRAITS_H
