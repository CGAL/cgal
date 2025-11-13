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

template <typename CircularKernel>
class Non_x_monotonic_Circular_arc_2 :
    public CircularKernel::Circular_arc_2 {
  using FT = typename CircularKernel::FT;
  using Point_2 = typename CircularKernel::Point_2;
  using Line_2 = typename CircularKernel::Line_2;
  using Circle_2 = typename CircularKernel::Circle_2;
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;

  using Base = typename CircularKernel::Circular_arc_2;

public:
  Non_x_monotonic_Circular_arc_2(): Base(){}

  Non_x_monotonic_Circular_arc_2(const Circle_2& c): Base(c){}
  // Not Documented
  Non_x_monotonic_Circular_arc_2(const Circle_2& support,
                                 const Line_2& l1, const bool b_l1,
                                 const Line_2& l2, const bool b_l2) :
    Base(support,l1,b_l1,l2,b_l2){}

  // Not Documented
  Non_x_monotonic_Circular_arc_2(const Circle_2& c,
                                 const Circle_2& c1, const bool b_1,
                                 const Circle_2& c2, const bool b_2) :
    Base(c,c1,b_1,c2,b_2)
  {}

  Non_x_monotonic_Circular_arc_2(const Point_2& start,
                                 const Point_2& middle,
                                 const Point_2& end) :
    Base(start,middle,end)
  {}

  Non_x_monotonic_Circular_arc_2(const Circle_2& support,
                                 const Circular_arc_point_2& begin,
                                 const Circular_arc_point_2& end) :
    Base(support,begin,end)
  {}

  Non_x_monotonic_Circular_arc_2(const Point_2& start,
                                 const Point_2& end,
                                 const FT& bulge) :
    Base(start,end,bulge)
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
  using Kernel = CircularKernel;
  using Curve_2 = internal::Non_x_monotonic_Circular_arc_2<CircularKernel>;
  using X_monotone_curve_2 = typename CircularKernel::Circular_arc_2;

  using Point = typename CircularKernel::Circular_arc_point_2;
  using Point_2 = typename CircularKernel::Circular_arc_point_2;

  using Multiplicity = std::size_t;

  using Has_left_category = CGAL::Tag_false;
  using Has_merge_category = CGAL::Tag_false;
  using Has_do_intersect_category = CGAL::Tag_false;

  using Left_side_category = Arr_oblivious_side_tag;
  using Bottom_side_category = Arr_oblivious_side_tag;
  using Top_side_category = Arr_oblivious_side_tag;
  using Right_side_category = Arr_oblivious_side_tag;

  Arr_circular_arc_traits_2(const CircularKernel& k = CircularKernel()) : ck(k) {}

  using Compare_x_2 = typename CircularKernel::Compare_x_2;
  using Compare_xy_2 = typename CircularKernel::Compare_xy_2;
  using Compare_y_at_x_2 = typename CircularKernel::Compare_y_at_x_2;
  using Compare_y_at_x_right_2 = typename CircularKernel::Compare_y_to_right_2;
  using Construct_max_vertex_2 = typename CircularKernel::Construct_circular_max_vertex_2;
  using Construct_min_vertex_2 = typename CircularKernel::Construct_circular_min_vertex_2;
  using Equal_2 = typename CircularKernel::Equal_2;
  // typedef typename CircularKernel::Make_x_monotone_2    Make_x_monotone_2;
  using Split_2 = typename CircularKernel::Split_2;
  using Intersect_2 = typename CircularKernel::Intersect_2;
  using Is_vertical_2 = typename CircularKernel::Is_vertical_2;

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
    OutputIterator operator()(const Curve_2& arc, OutputIterator oi) const {
      using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;

      std::vector<Make_x_monotone_result> objs;
      CircularKernel().make_x_monotone_2_object()(arc, std::back_inserter(objs));
      for (const auto& obj : objs) {
        if (const auto* p = std::get_if<Point_2>(&obj)) {
          *oi++ = Make_x_monotone_result(*p);
          continue;
        }
        if (const auto* xcv = std::get_if<X_monotone_curve_2>(&obj)) {
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

#endif
