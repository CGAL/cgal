// Copyright (c) 2011 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri


#ifndef CGAL_INTERNAL_STATIC_FILTERS_COMPARE_Y_2_H
#define CGAL_INTERNAL_STATIC_FILTERS_COMPARE_Y_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/tools.h>
#include <cmath>
#include <iostream>

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {


template < typename K_base >
class Compare_y_2
  : public K_base::Compare_y_2
{
  typedef typename K_base::Point_2   Point_2;
  typedef typename K_base::Line_2    Line_2;
  typedef typename K_base::Compare_y_2   Base;

public:

  typedef typename Base::result_type  result_type;

  using Base::operator();

  result_type operator()(const Point_2 &p, const Point_2& q) const
  {
    CGAL_BRANCH_PROFILER(std::string("semi-static attempts/calls to   : ") +
                         std::string(CGAL_PRETTY_FUNCTION), tmp);

    Get_approx<Point_2> get_approx; // Identity functor for all points
                                    // but lazy points
    double py, qy;

    if (fit_in_double(get_approx(p).y(), py) && fit_in_double(get_approx(q).y(), qy) )
    {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      return compare(py, qy);
    }
    return Base::operator()(p, q);
  }



}; // end class Compare_y_2

} // end namespace Static_filters_predicates

} // end namespace internal

} // end namespace CGAL

#endif  // CGAL_INTERNAL_STATIC_FILTERS_COMPARE_Y_2_H
