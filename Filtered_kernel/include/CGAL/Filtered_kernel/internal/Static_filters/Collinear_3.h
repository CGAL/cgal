// Copyright (c) 2017 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_INTERNAL_STATIC_FILTERS_COLLINEAR_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_COLLINEAR_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Static_filter_error.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/tools.h>
#include <cmath>

namespace CGAL { namespace internal { namespace Static_filters_predicates {

template < typename K_base >
class Collinear_3
  : public K_base::Collinear_3
{
  typedef typename K_base::Point_3          Point_3;
  typedef typename K_base::Vector_3         Vector_3;
  typedef typename K_base::Sphere_3         Sphere_3;
  typedef typename K_base::Tetrahedron_3    Tetrahedron_3;
  typedef typename K_base::Collinear_3    Base;

public:

  typedef typename Base::result_type  result_type;
  result_type
  operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r) const
  {
    CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Collinear_3", tmp);

    Get_approx<Point_3> get_approx; // Identity functor for all points
    // but lazy points.

    double px, py, pz, qx, qy, qz, rx, ry, rz;

    if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
        fit_in_double(get_approx(p).z(), pz) &&
        fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
        fit_in_double(get_approx(q).z(), qz) &&
        fit_in_double(get_approx(r).x(), rx) && fit_in_double(get_approx(r).y(), ry) &&
        fit_in_double(get_approx(r).z(), rz))
    {
      CGAL_BRANCH_PROFILER_BRANCH_1(tmp);

      double dpx = (px - rx);
      double dqx = (qx - rx);
      double dpy = (py - ry);
      double dqy = (qy - ry);

      double double_tmp_result = ((dpx * dqy) - (dpy * dqx));

      double max1 = CGAL::abs(dpx);
      if (max1 < CGAL::abs(dqx))
        max1 = CGAL::abs(dqx);

      double max2 = CGAL::abs(dpy);
      if (max2 < CGAL::abs(dqy))
        max2 = CGAL::abs(dqy);

      double lower_bound_1 = max2;
      double upper_bound_1 = max2;
      if (max1 < lower_bound_1)
        lower_bound_1 = max1;
      else if (max1 > upper_bound_1)
        upper_bound_1 = max1;

      if ((lower_bound_1 >= 5.00368081960964635413e-147) && (upper_bound_1 <= 1.67597599124282407923e+153)) {
        double eps = (8.88720573725927976811e-16 * (max1 * max2));
        if (CGAL::abs(double_tmp_result) > eps )
          return false;
      }

      double dpz = (pz - rz);
      double dqz = (qz - rz);

      double double_tmp_result_k3Lzf6g = ((dpx * dqz) - (dpz * dqx));

      double max3 = CGAL::abs(dpz);
      if (max3 < CGAL::abs(dqz))
        max3 = CGAL::abs(dqz);

      lower_bound_1 = max1;
      upper_bound_1 = max1;

      if (max3 < lower_bound_1)
        lower_bound_1 = max3;
      else if (max3 > upper_bound_1)
        upper_bound_1 = max3;

      if ((lower_bound_1 >= 5.00368081960964635413e-147) && (upper_bound_1 <= 1.67597599124282407923e+153)) {
        double eps = (8.88720573725927976811e-16 * (max1 * max3));
        if (CGAL::abs(double_tmp_result_k3Lzf6g) > eps)
          return false;
      }

      double double_tmp_result_AvrrXBP = ((dpy * dqz) - (dpz * dqy));

      lower_bound_1 = max2;
      upper_bound_1 = max2;

      if (max3 < lower_bound_1)
        lower_bound_1 = max3;
      else if (max3 > upper_bound_1)
        upper_bound_1 = max3;

      if ((lower_bound_1 >= 5.00368081960964635413e-147) && (upper_bound_1 <= 1.67597599124282407923e+153)) {
        double eps = (8.88720573725927976811e-16 * (max2 * max3));
        if (CGAL::abs(double_tmp_result_AvrrXBP) > eps)
          return false;
      }

      CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
    }

    return Base::operator()(p, q, r);
  }

};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_COLLINEAR_3_H
