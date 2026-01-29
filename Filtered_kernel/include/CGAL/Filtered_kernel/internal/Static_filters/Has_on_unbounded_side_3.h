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


#ifndef CGAL_INTERNAL_STATIC_FILTERS_HAS_ON_UNBOUNDED_SIDE_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_HAS_ON_UNBOUNDED_SIDE_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/tools.h>
#include <cmath>
#include <iostream>

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {


template < typename K_base >
class Has_on_unbounded_side_3
  : public K_base::Has_on_unbounded_side_3
{
  typedef typename K_base::Boolean                 Boolean;
  typedef typename K_base::Point_3                 Point_3;
  typedef typename K_base::Sphere_3                Sphere_3;
  typedef typename K_base::Iso_cuboid_3            Iso_cuboid_3;
  typedef typename K_base::Has_on_unbounded_side_3 Base;

public:
  using Base::operator();

  Boolean
  operator()(const Sphere_3& s, const Iso_cuboid_3& b) const
  {
    CGAL_BRANCH_PROFILER_3(std::string("semi-static failures/attempts/calls to   : ") +
      std::string(CGAL_PRETTY_FUNCTION), tmp);

    Get_approx<Point_3> get_approx; // Identity functor for all points
    const Point_3& c = s.center();

    double scx, scy, scz, ssr;
    double bxmin = b.xmin(), bymin = b.ymin(), bzmin = b.zmin(),
      bxmax = b.xmax(), bymax = b.ymax(), bzmax = b.zmax();

    if (fit_in_double(get_approx(c).x(), scx) &&
      fit_in_double(get_approx(c).y(), scy) &&
      fit_in_double(get_approx(c).z(), scz) &&
      fit_in_double(s.squared_radius(), ssr))
    {
      CGAL_BRANCH_PROFILER_BRANCH_1(tmp);

      if ((ssr < 1.11261183279326254436e-293) || (ssr > 2.80889552322236673473e+306)) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(s, b);
      }
      double distance = 0;
      double max1 = 0;
      double double_tmp_result = 0;
      double eps = 0;
      if (scx < bxmin)
      {
        double bxmin_scx = bxmin - scx;
        max1 = bxmin_scx;

        distance = square(bxmin_scx);
        double_tmp_result = (distance - ssr);

        if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153)) {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(s, b);
        }

        eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

        if (double_tmp_result > eps) {
          return true;
        }
      }
      else if (scx > bxmax)
      {
        double scx_bxmax = scx - bxmax;
        max1 = scx_bxmax;

        distance = square(scx_bxmax);
        double_tmp_result = (distance - ssr);

        if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153)) {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(s, b);
        }

        eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

        if (double_tmp_result > eps) {
          return true;
        }
      }
      else
      {
        double scx_bxmin = scx - bxmin;
        double bxmax_scx = bxmax - scx;
        max1 = (std::min)(scx_bxmin, bxmax_scx);

        distance = square(max1);
        double_tmp_result = (distance - ssr);

        if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153)) {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(s, b);
        }
      }


      if (scy < bymin)
      {
        double bymin_scy = bymin - scy;
        if (max1 < bymin_scy) {
          max1 = bymin_scy;
        }

        distance += square(bymin_scy);
      }
      else if (scy > bymax)
      {
        double scy_bymax = scy - bymax;
        if (max1 < scy_bymax) {
          max1 = scy_bymax;
        }
        distance += square(scy_bymax);
      }
      else {
        double scy_bymin = scy - bymin;
        double bymax_scy = bymax - scy;
        double d = (std::min)(scy_bymin, bymax_scy);
        if (max1 < d) {
          max1 = d;
        }

        distance += square(d);
      }

      double_tmp_result = (distance - ssr);

      if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153)) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(s, b);
      }

      eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

      if (double_tmp_result > eps) {
        return true;
      }


      if (scz < bzmin)
      {
        double bzmin_scz = bzmin - scz;
        if (max1 < bzmin_scz) {
          max1 = bzmin_scz;
        }
        distance += square(bzmin_scz);
      }
      else if (scz > bzmax)
      {
        double scz_bzmax = scz - bzmax;
        if (max1 < scz_bzmax) {
          max1 = scz_bzmax;
        }

        distance += square(scz_bzmax);
      }
      else
      {
        double scz_bzmin = scz - bzmin;
        double bzmax_scz = bzmax - scz;
        double d = (std::min)(scz_bzmin, bzmax_scz);

        if (max1 < d) {
          max1 = d;
        }

        distance += square(d);
      }

      double_tmp_result = (distance - ssr);
      if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153)) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(s, b);
      }

      eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

      if (double_tmp_result > eps) {
        return true;
      }
      else if (double_tmp_result <= -eps) {
        return false;
      }
      else {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(s, b);
      }

      CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
    }
    return Base::operator()(s, b);
  }

  Boolean
  operator()(const Iso_cuboid_3& b, const Sphere_3& s) const
  {
    this->operator()(s, b);
  }


}; // end class Has_on_unbounded_side_3

} // end namespace Static_filters_predicates

} // end namespace internal

} // end namespace CGAL

#endif  // CGAL_INTERNAL_STATIC_FILTERS_HAS_ON_UNBOUNDED_SIDE_3_H
