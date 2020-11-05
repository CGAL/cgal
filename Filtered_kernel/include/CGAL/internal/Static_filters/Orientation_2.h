// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_2_H
#define CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_2_H

#include <CGAL/Profile_counter.h>
#include <CGAL/determinant.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>

#include <cmath>

namespace CGAL { namespace internal { namespace Static_filters_predicates {

template < typename K_base >
class Orientation_2
  : public K_base::Orientation_2
{
  typedef typename K_base::Point_2          Point_2;
  typedef typename K_base::Vector_2         Vector_2;
  typedef typename K_base::Circle_2          Circle_2;

  typedef typename K_base::Orientation_2    Base;

public:

  typedef typename Base::result_type  result_type;

  using Base::operator();

  Orientation
  operator()(const Point_2 &p, const Point_2 &q, const Point_2 &r) const

  {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Orientation_2", tmp);

      double px, py, qx, qy, rx, ry;

      if (fit_in_double(p.x(), px) && fit_in_double(p.y(), py) &&
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy) &&
          fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry))
      {
          CGAL_BRANCH_PROFILER_BRANCH_1(tmp);

          double pqx = qx - px;
          double pqy = qy - py;
          double prx = rx - px;
          double pry = ry - py;

          double det = CGAL::determinant(pqx, pqy,
                                         prx, pry);

          // Then semi-static filter.
          double maxx = CGAL::abs(pqx);
          double maxy = CGAL::abs(pqy);

          double aprx = CGAL::abs(prx);
          double apry = CGAL::abs(pry);

          if (maxx < aprx) maxx = aprx;

          if (maxy < apry) maxy = apry;

          // Sort them
          if (maxx > maxy)  std::swap(maxx, maxy);

          // Protect against underflow in the computation of eps.
          if (maxx < 1e-146) /* sqrt(min_double/eps) */ {
            if (maxx == 0)
              return ZERO;
          }
          // Protect against overflow in the computation of det.
          else if (maxy < 1e153) /* sqrt(max_double [hadamard]/2) */ {
            double eps = 8.8872057372592798e-16 * maxx * maxy;
            if (det > eps)  return POSITIVE;
            if (det < -eps) return NEGATIVE;
          }

          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
      }

      return Base::operator()(p, q, r);
  }

  // Computes the epsilon for Orientation_2.
  static double compute_epsilon()
  {
    typedef Static_filter_error F;
    F t1 = F(1, F::ulp()/2);         // First translation
    F det = CGAL::determinant(t1, t1,
                              t1, t1); // Full det
    double err = det.error();
    err += err * 2 * F::ulp(); // Correction due to "epsilon * maxx * maxy".
    std::cerr << "*** epsilon for Orientation_2 = " << err << std::endl;
    return err;
  }
};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_2_H
