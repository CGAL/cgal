// Copyright (c) 2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Stéphane Tayeb, Pierre Alliez, Camille Wormser
//

#ifndef CGAL_AABB_STATICALLY_FILTERED_TRAITS_3_H
#define CGAL_AABB_STATICALLY_FILTERED_TRAITS_3_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/Filtered_kernel/internal/Static_filters/Static_filter_error.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/tools.h>
#include <CGAL/AABB_tree/internal/AABB_filtered_traits_3.h>

namespace CGAL {

template<typename GeomTraits, typename AABBPrimitive, typename BboxMap = Default>
class AABB_statically_filtered_traits_3 : public AABB_filtered_traits_base_3<GeomTraits, AABBPrimitive, BboxMap> {
  using Base = AABB_filtered_traits_base_3<GeomTraits, AABBPrimitive, BboxMap>;

public:
  class Compare_distance : public Base::Compare_distance {
  public:
    using Point = typename GeomTraits::Point_3;
    using Bounding_box = Bbox_3;
    using Sphere_3 = typename GeomTraits::Sphere_3;

    Comparison_result operator()(const Point& p, const Bounding_box& b, const Point& bound) const {
      Sphere_3 s = GeomTraits().construct_sphere_3_object()(p, GeomTraits().compute_squared_distance_3_object()(p, bound));
      CGAL_BRANCH_PROFILER_3(std::string("semi-static failures/attempts/calls to   : ") +
        std::string(CGAL_PRETTY_FUNCTION), tmp);

      internal::Static_filters_predicates::Get_approx<Point> get_approx; // Identity functor for all points
      const Point& c = s.center();

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

          if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153))
            return CGAL::SMALLER;

          eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

          if (double_tmp_result > eps)
            return CGAL::LARGER;
        }
        else if (scx > bxmax)
        {
          double scx_bxmax = scx - bxmax;
          max1 = scx_bxmax;

          distance = square(scx_bxmax);
          double_tmp_result = (distance - ssr);

          if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153))
            return CGAL::SMALLER;

          eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

          if (double_tmp_result > eps)
            return CGAL::LARGER;
        }


        if (scy < bymin)
        {
          double bymin_scy = bymin - scy;
          if (max1 < bymin_scy) {
            max1 = bymin_scy;
          }

          distance += square(bymin_scy);
          double_tmp_result = (distance - ssr);

          if ((max1 < 3.33558365626356687717e-147) || ((max1 > 1.67597599124282407923e+153)))
            return CGAL::SMALLER;

          eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

          if (double_tmp_result > eps)
            return CGAL::LARGER;
        }
        else if (scy > bymax)
        {
          double scy_bymax = scy - bymax;
          if (max1 < scy_bymax) {
            max1 = scy_bymax;
          }
          distance += square(scy_bymax);
          double_tmp_result = (distance - ssr);

          if (((max1 < 3.33558365626356687717e-147)) || ((max1 > 1.67597599124282407923e+153)))
            return CGAL::SMALLER;

          eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

          if (double_tmp_result > eps)
            return CGAL::LARGER;
        }


        if (scz < bzmin)
        {
          double bzmin_scz = bzmin - scz;
          if (max1 < bzmin_scz) {
            max1 = bzmin_scz;
          }
          distance += square(bzmin_scz);
          double_tmp_result = (distance - ssr);

          if (((max1 < 3.33558365626356687717e-147)) || ((max1 > 1.67597599124282407923e+153)))
            return CGAL::SMALLER;

          eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

          if (double_tmp_result > eps)
            return CGAL::LARGER;
        }
        else if (scz > bzmax)
        {
          double scz_bzmax = scz - bzmax;
          if (max1 < scz_bzmax) {
            max1 = scz_bzmax;
          }

          distance += square(scz_bzmax);
          double_tmp_result = (distance - ssr);

          if (((max1 < 3.33558365626356687717e-147)) || ((max1 > 1.67597599124282407923e+153)))
            return CGAL::SMALLER;

          eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

          if (double_tmp_result > eps)
            return CGAL::LARGER;
        }

        // double_tmp_result and eps were growing all the time
        // no need to test for > eps as done earlier in at least one case
        if (double_tmp_result < -eps)
          return CGAL::LARGER;
        else
          return CGAL::SMALLER;

        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
      }
      return Base::operator()(s, b);
    }
  };
  AABB_statically_filtered_traits_3() : Base() {}
  AABB_statically_filtered_traits_3(BboxMap bbm) : Base(bbm) {}
};

} //namespace CGAL

#endif
