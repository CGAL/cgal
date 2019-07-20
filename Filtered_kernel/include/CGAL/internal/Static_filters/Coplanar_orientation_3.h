// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_INTERNAL_STATIC_FILTERS_COPLANAR_ORIENTATION_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_COPLANAR_ORIENTATION_3_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Profile_counter.h>

namespace CGAL { namespace internal { namespace Static_filters_predicates {

// XXX :
// Is this really useful to optimize this predicate ?
// I de-prioritize it for now.


// This one is easy : it's just 3 Orientation_2.

template < typename P3 >
struct Point_23_xy {
    const P3 &p;

    Point_23_xy(const P3& pp) : p(pp) {}

    const FT& x() const { return p.x(); }
    const FT& y() const { return p.y(); }
};

template < typename Kernel >
class Coplanar_orientation_3
  : public Kernel::Coplanar_orientation_3
{
  typedef typename Kernel::Point_3                  Point_3;
  typedef typename Kernel::Coplanar_orientation_3   Base;

public:
  typedef Orientation result_type;

  Orientation operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r) const
  {
      return opti_coplanar_orientationC3(
	    to_double(p.x()), to_double(p.y()), to_double(p.z()),
	    to_double(q.x()), to_double(q.y()), to_double(q.z()),
	    to_double(r.x()), to_double(r.y()), to_double(r.z()));
  }

  Orientation operator()(const Point_3 &p, const Point_3 &q,
                         const Point_3 &r, const Point_3 &s) const
  {
      return opti_coplanar_orientationC3(
	    to_double(p.x()), to_double(p.y()), to_double(p.z()),
	    to_double(q.x()), to_double(q.y()), to_double(q.z()),
	    to_double(r.x()), to_double(r.y()), to_double(r.z()),
	    to_double(s.x()), to_double(s.y()), to_double(s.z()));
  }

private:
  Orientation
  opti_coplanar_orientationC3(double px, double py, double pz,
                              double qx, double qy, double qz,
		              double rx, double ry, double rz) const
  {
      CGAL_PROFILER("Coplanar_orientation_3 #1 calls");

      std::pair<bool, Orientation> oxy_pqr = orient_2d(px,py,qx,qy,rx,ry);
      if (oxy_pqr != COLLINEAR)
          return oxy_pqr;

      CGAL_PROFILER("Coplanar_orientation_3 #1 step2");

      Orientation oyz_pqr = orient_2d(py,pz,qy,qz,ry,rz);
      if (oyz_pqr != COLLINEAR)
          return oyz_pqr;

      CGAL_PROFILER("Coplanar_orientation_3 #1 step3");

      return orient_2d(px,pz,qx,qz,rx,rz);

  }

  Orientation
  opti_coplanar_orientationC3(double px, double py, double pz,
                              double qx, double qy, double qz,
		              double rx, double ry, double rz,
		              double sx, double sy, double sz) const
  {
      CGAL_PROFILER("Coplanar_orientation_3 #2 calls");

      Orientation oxy_pqr = orient_2d(px,py,qx,qy,rx,ry);
      if (oxy_pqr != COLLINEAR)
          return Orientation( oxy_pqr *
		              orient_2d(px,py,qx,qy,sx,sy));

      CGAL_PROFILER("Coplanar_orientation_3 #2 step2");

      Orientation oyz_pqr = orient_2d(py,pz,qy,qz,ry,rz);
      if (oyz_pqr != COLLINEAR)
          return Orientation( oyz_pqr *
		              orient_2d(py,pz,qy,qz,sy,sz));

      CGAL_PROFILER("Coplanar_orientation_3 #2 step3");

      Orientation oxz_pqr = orient_2d(px,pz,qx,qz,rx,rz);
      CGAL_kernel_assertion(oxz_pqr != COLLINEAR);
      return Orientation( oxz_pqr * orient_2d(px,pz,qx,qz,sx,sz));
  }

  // FIXME : Some code duplicated from Orientation_2...
  std::pair<bool, Orientation>
  orient_2d(double px, double py, double qx, double qy, double rx, double ry) const
  {
    /*
    double px = p.x();
    double py = p.y();
    double qx = q.x();
    double qy = q.y();
    double rx = r.x();
    double ry = r.y();
    */

    CGAL_PROFILER("orient2d calls");

    double pqx = qx-px;
    double pqy = qy-py;
    double prx = rx-px;
    double pry = ry-py;

    double det = determinant(pqx, pqy,
                                   prx, pry);

    // Then semi-static filter.
    double maxx = CGAL::abs(px);
    double maxy = CGAL::abs(py);

    double aqx = CGAL::abs(qx);
    double aqy = CGAL::abs(qy);

    double arx = CGAL::abs(rx);
    double ary = CGAL::abs(ry);

    if (maxx < aqx) maxx = aqx;
    if (maxx < arx) maxx = arx;
    if (maxy < aqy) maxy = aqy;
    if (maxy < ary) maxy = ary;
    double eps = 3.55271e-15 * maxx * maxy;

    if (det > eps)  return std::make_pair(true, POSITIVE);
    if (det < -eps) return std::make_pair(true, NEGATIVE);

    CGAL_PROFILER("orient2d semi-static failures");

    return std::make_pair(false, ZERO);
  }

};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_COPLANAR_ORIENTATION_3_H
