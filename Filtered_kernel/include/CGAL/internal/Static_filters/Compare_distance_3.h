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


#ifndef CGAL_INTERNAL_STATIC_FILTERS_COMPARE_DISTANCE_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_COMPARE_DISTANCE_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/tools.h>
#include <cmath>
#include <iostream>

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {


template < typename K_base >
class Compare_distance_3
  : public K_base::Compare_distance_3
{
  typedef typename K_base::Point_3   Point_3;
  typedef typename K_base::Vector_3  Vector_3;
  typedef typename K_base::Compare_distance_3   Base;

public:

  typedef typename Base::result_type  result_type;

  using Base::operator();

  result_type operator()(const Point_3 &p, const Point_3& q, const Point_3& r) const
  {
    CGAL_BRANCH_PROFILER(std::string("semi-static attempts/calls to   : ") +
                         std::string(CGAL_PRETTY_FUNCTION), tmp);

    Get_approx<Point_3> get_approx; // Identity functor for all points
                                    // but lazy points

    if(q == r){
      return EQUAL;
    }
    double px, py, pz, qx, qy, qz, rx, ry, rz;

    if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
        fit_in_double(get_approx(p).z(), pz) &&
        fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
        fit_in_double(get_approx(q).z(), qz) &&
        fit_in_double(get_approx(r).x(), rx) && fit_in_double(get_approx(r).y(), ry) &&
        fit_in_double(get_approx(r).z(), rz))
    {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      double qpx;
      qpx = (qx - px);
      double qpy;
      qpy = (qy - py);
      double qpz;
      qpz = (qz - pz);
      double qp2;
      qp2 = ((square( qpx ) + square( qpy )) + square( qpz ));
      double rpx;
      rpx = (rx - px);
      double rpy;
      rpy = (ry - py);
      double rpz;
      rpz = (rz - pz);
      double rp2;
      rp2 = ((square( rpx ) + square( rpy )) + square( rpz ));
      Sign int_tmp_result = EQUAL;
      double double_tmp_result;
      double eps;
      double_tmp_result = (qp2 - rp2);
      double max1 = CGAL::abs(qpx);
      if( (max1 < CGAL::abs(qpy)) )
        {
          max1 = CGAL::abs(qpy);
        }
      if( (max1 < CGAL::abs(qpz)) )
        {
          max1 = CGAL::abs(qpz);
        }
      if( (max1 < CGAL::abs(rpx)) )
        {
          max1 = CGAL::abs(rpx);
        }
      if( (max1 < CGAL::abs(rpy)) )
        {
          max1 = CGAL::abs(rpy);
        }
      if( (max1 < CGAL::abs(rpz)) )
        {
          max1 = CGAL::abs(rpz);
        }
      if( (max1 < 2.42701102401884262260e-147) )
        {
          return Base::operator()(p, q, r);
        }
      else
        {
          if( (max1 > 8.37987995621411946582e+152) )
            {
              return Base::operator()(p, q, r);
            }
          eps = (3.77746921267322435884e-15 * (max1 * max1));
          if( (double_tmp_result > eps) )
            {
              int_tmp_result = LARGER;
            }
          else
            {
              if( (double_tmp_result < -eps) )
                {
                  int_tmp_result = SMALLER;
                }
              else
                {
                  return Base::operator()(p, q, r);
                }
            }
        }
      return int_tmp_result;
    }
    return Base::operator()(p, q, r);
  }


}; // end class Compare_distance_3

} // end namespace Static_filters_predicates

} // end namespace internal

} // end namespace CGAL

#endif  // CGAL_INTERNAL_STATIC_FILTERS_COMPARE_DISTANCE_3_H
