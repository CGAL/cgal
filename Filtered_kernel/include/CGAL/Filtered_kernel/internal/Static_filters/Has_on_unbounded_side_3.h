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
      bool center_inside =  ! ( (scx < bxmin) || (scx > bxmax) || (scy < bymin) || (scy > bymax)|| (scz < bzmin) || (scz > bzmax) );
      if(center_inside){
        return false;
      }

      //  AF: What does this guard against ?
      if ((ssr < 1.11261183279326254436e-293) || (ssr > 2.80889552322236673473e+306)) {
        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
        return Base::operator()(s, b);
      }

      // When the center is outside, the closest point has the following coordinates
      double px = std::clamp(scx, bxmin, bxmax); // equivalent to:  std::min(std::max(scx,bxmin),bxmax);
      double py = std::clamp(scy, bymin, bymax);
      double pz = std::clamp(scz, bzmin, bzmax);

      // AF:  We now need an error bound on the folling expression
      double sqdist = CGAL::square(px-scx) + CGAL::square(py-scy) + CGAL::square(pz-scz);
      double eps = 0;
      double delta = sqdist - ssr;
      if(delta > eps || delta < -eps){
        return ssr < sqdist;
      }

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
