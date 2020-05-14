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


#ifndef CGAL_INTERNAL_STATIC_FILTERS_EQUAL_2_H
#define CGAL_INTERNAL_STATIC_FILTERS_EQUAL_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/tools.h>
#include <cmath>
#include <iostream>

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {


template < typename K_base >
class Equal_2
  : public K_base::Equal_2
{
  typedef typename K_base::FT        FT;
  typedef typename K_base::Point_2   Point_2;
  typedef typename K_base::Vector_2  Vector_2;
  typedef typename K_base::Equal_2   Base;

public:

  typedef typename Base::result_type  result_type;


#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else // CGAL_CFG_MATCHING_BUG_6
  template <typename T>
  result_type
  operator()(const T& t1, const T& t2) const
  {
    return Base()(t1,t2);
  }
#endif // CGAL_CFG_MATCHING_BUG_6


  result_type operator()(const Point_2 &p, const Point_2& q) const
  {
    CGAL_BRANCH_PROFILER(std::string("semi-static attempts/calls to   : ") +
                         std::string(CGAL_PRETTY_FUNCTION), tmp);

    double px, py, qx, qy;
    init_double(px, py,qx, qy, (FT*)(0));

    if (fit_in_double(p.x(), px) && fit_in_double(p.y(), py) &&
        fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy) )
    {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      return px == qx && py == qy;
    }
    return Base::operator()(p, q);
  }


  result_type operator()(const Vector_2 &p, const Vector_2& q) const
  {
    CGAL_BRANCH_PROFILER(std::string("semi-static attempts/calls to   : ") +
                         std::string(CGAL_PRETTY_FUNCTION), tmp);

    double px, py, qx, qy;

    if (fit_in_double(p.x(), px) && fit_in_double(p.y(), py) &&
        fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy) )
    {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      return px == qx && py == qy;
    }
    return Base::operator()(p, q);
  }

}; // end class Equal_2

} // end namespace Static_filters_predicates

} // end namespace internal

} // end namespace CGAL

#endif  // CGAL_INTERNAL_STATIC_FILTERS_EQUAL_3_H
