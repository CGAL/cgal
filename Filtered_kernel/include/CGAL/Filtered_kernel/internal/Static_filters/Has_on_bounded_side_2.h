// Copyright (c) 2026 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau


#ifndef CGAL_INTERNAL_STATIC_FILTERS_HAS_ON_BOUNDED_SIDE_2_H
#define CGAL_INTERNAL_STATIC_FILTERS_HAS_ON_BOUNDED_SIDE_2_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/tools.h>
#include <cmath>
#include <iostream>

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {


template < typename K_base >
class Has_on_bounded_side_2
  : public K_base::Has_on_bounded_side_2
{
  typedef typename K_base::Boolean               Boolean;
  typedef typename K_base::Point_2               Point_2;
  typedef typename K_base::Circle_2              Circle_2;
  typedef typename K_base::Segment_2             Segment_2;
  typedef typename K_base::Has_on_bounded_side_2 Base;

public:
  using Base::operator();

  Boolean
  operator()(const Circle_2& c, const Segment_2& s) const
  {
    CGAL_BRANCH_PROFILER_3(std::string("semi-static failures/attempts/calls to   : ") +
      std::string(CGAL_PRETTY_FUNCTION), tmp);

    Get_approx<typename K_base::Point_2> get_approx; // Identity functor for all points
    const Point_2& center = c.center();
    const Point_2& src = s.source();
    const Point_2& tgt = s.target();

    double cx, cy, csr;
    double sx, sy, tx, ty;

    if (fit_in_double(get_approx(center).x(), cx) &&
      fit_in_double(get_approx(center).y(), cy) &&
      fit_in_double(get_approx(src).x(), sx) &&
      fit_in_double(get_approx(src).y(), sy) &&
      fit_in_double(get_approx(tgt).x(), tx) &&
      fit_in_double(get_approx(tgt).y(), ty) &&
      fit_in_double(c.squared_radius(), csr))
    {
      CGAL_BRANCH_PROFILER_BRANCH_1(tmp);

      if ((csr < 1.11261183279326254436e-293) || (csr > 2.80889552322236673473e+306)) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(c, s);
      }

      double distance = 0;
      double max1 = 0;
      double double_tmp_result = 0;
      double eps = 0;

      double dx = cx - sx;
      max1 = dx;
      distance = square(dx);
      double_tmp_result = (distance - csr);

      if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153)) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(c, s);
      }

      eps = 1.99986535548615598560e-15 * (std::max)(csr, square(max1));

      if (double_tmp_result > eps) {
        return false;
      }

      double dy = cy - sy;
      if (max1 < dy) {
        max1 = dy;
      }
      distance += square(dy);
      double_tmp_result = (distance - csr);
      if ((max1 < 3.33558365626356687717e-147) || ((max1 > 1.67597599124282407923e+153))) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(c, s);
      }
      eps = 1.99986535548615598560e-15 * (std::max)(csr, square(max1));

      if (double_tmp_result > eps) {
        return false;
      }
      else if (double_tmp_result > -eps) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(c, s);
      }

      distance = max1 = double_tmp_result = eps = 0;

      dx = cx - tx;
      max1 = dx;
      distance = square(dx);
      double_tmp_result = (distance - csr);

      if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153)) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(c, s);
      }

      eps = 1.99986535548615598560e-15 * (std::max)(csr, square(max1));

      if (double_tmp_result > eps) {
        return false;
      }

      dy = cy - ty;
      if (max1 < dy) {
        max1 = dy;
      }
      distance += square(dy);
      double_tmp_result = (distance - csr);
      if ((max1 < 3.33558365626356687717e-147) || ((max1 > 1.67597599124282407923e+153))) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(c, s);
      }
      eps = 1.99986535548615598560e-15 * (std::max)(csr, square(max1));

      if (double_tmp_result > eps) {
        return false;
      }
      else if (double_tmp_result > -eps) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(c, s);
      }
      else return true;

      CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
    }
    return Base::operator()(c, s);
  }

  Boolean
  operator()(const Segment_2& s, const Circle_2& c) const
  {
    this->operator()(c, s);
  }


}; // end class Has_on_bounded_side_2

} // end namespace Static_filters_predicates

} // end namespace internal

} // end namespace CGAL

#endif  // CGAL_INTERNAL_STATIC_FILTERS_HAS_ON_BOUNDED_SIDE_2_H
