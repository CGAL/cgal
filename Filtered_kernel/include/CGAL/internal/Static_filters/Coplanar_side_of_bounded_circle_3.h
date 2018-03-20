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

#ifndef CGAL_INTERNAL_STATIC_FILTERS_COPLANAR_SIDE_OF_BOUNDED_CIRCLE_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_COPLANAR_SIDE_OF_BOUNDED_CIRCLE_3_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Profile_counter.h>
// #include <CGAL/internal/Static_filters/Static_filter_error.h> // Only used to precompute constants

namespace CGAL { namespace internal { namespace Static_filters_predicates {

template <class Point>
class Side_of_bounded_circle_3
{
  // Computes the epsilon for In_circle_3.
  static double cir_3()
  {
      // NOTE : This produces a buggy degree warning.
      // Maybe there's a bug in the formula.
    typedef CGAL::Static_filter_error F;
    F t1 = F(1)-F(1);         // First translation
    F sq = t1*t1+t1*t1+t1*t1; // squares
    F n1 = t1*t1 - t1*t1;     // normal vector
    F sq_n1 = n1*n1 + n1*n1 + n1*n1;
    F det = determinant(t1, t1, t1, sq,
                              t1, t1, t1, sq,
                              t1, t1, t1, sq,
                              n1, n1, n1, sq_n1); // Full det
    double err = det.error();
    std::cerr << "*** epsilon for In_circle_3 = " << err << std::endl;
    return err;
  }

public:
  typedef Bounded_side result_type;

  Bounded_side operator()(const Point &p, const Point &q,
	                  const Point &r, const Point &t) const
  {
    return opti_coplanar_in_circleC3(
	    to_double(p.x()), to_double(p.y()), to_double(p.z()),
            to_double(q.x()), to_double(q.y()), to_double(q.z()),
	    to_double(r.x()), to_double(r.y()), to_double(r.z()),
	    to_double(t.x()), to_double(t.y()), to_double(t.z()));
  }

  Bounded_side
  opti_coplanar_in_circleC3(double px, double py, double pz,
                            double qx, double qy, double qz,
		            double rx, double ry, double rz,
		            double tx, double ty, double tz) const
  {
    CGAL_PROFILER("In_circle_3 calls");

    double ptx = px - tx;
    double pty = py - ty;
    double ptz = pz - tz;
    double pt2 = CGAL_NTS square(ptx) + CGAL_NTS square(pty) +
	         CGAL_NTS square(ptz);
    double qtx = qx - tx;
    double qty = qy - ty;
    double qtz = qz - tz;
    double qt2 = CGAL_NTS square(qtx) + CGAL_NTS square(qty) +
	         CGAL_NTS square(qtz);
    double rtx = rx - tx;
    double rty = ry - ty;
    double rtz = rz - tz;
    double rt2 = CGAL_NTS square(rtx) + CGAL_NTS square(rty) +
	         CGAL_NTS square(rtz);
    double pqx = qx - px;
    double pqy = qy - py;
    double pqz = qz - pz;
    double prx = rx - px;
    double pry = ry - py;
    double prz = rz - pz;
    double vx = pqy*prz - pqz*pry;
    double vy = pqz*prx - pqx*prz;
    double vz = pqx*pry - pqy*prx;
    double v2 = CGAL_NTS square(vx) + CGAL_NTS square(vy) +
	        CGAL_NTS square(vz);

    double det = determinant(ptx,pty,ptz,pt2,
                                   rtx,rty,rtz,rt2,
                                   qtx,qty,qtz,qt2,
                                   vx,vy,vz,v2);

    // Compute the semi-static bound.
    double maxx = CGAL::abs(px);
    double maxy = CGAL::abs(py);
    double maxz = CGAL::abs(pz);

    double aqx = CGAL::abs(qx);
    double aqy = CGAL::abs(qy);
    double aqz = CGAL::abs(qz);

    double arx = CGAL::abs(rx);
    double ary = CGAL::abs(ry);
    double arz = CGAL::abs(rz);


    double atx = CGAL::abs(tx);
    double aty = CGAL::abs(ty);
    double atz = CGAL::abs(tz);

    if (maxx < aqx) maxx = aqx;
    if (maxx < arx) maxx = arx;
    if (maxx < atx) maxx = atx;

    if (maxy < aqy) maxy = aqy;
    if (maxy < ary) maxy = ary;
    if (maxy < aty) maxy = aty;

    if (maxz < aqz) maxz = aqz;
    if (maxz < arz) maxz = arz;
    if (maxz < atz) maxz = atz;

    double d = std::max(maxx, std::max(maxy, maxz));
    double eps = 3.27418e-11 * d * d * d * d * d * d;

    if (det >  eps) return ON_BOUNDED_SIDE;
    if (det < -eps) return ON_UNBOUNDED_SIDE;

    CGAL_PROFILER("In_circle_3 semi-static failures");

    typedef Filtered_kernel<Simple_cartesian<double>, MP_Float> K;
    typedef K::Point_3 P;

    return coplanar_side_of_bounded_circle(P(px,py,pz), P(qx,qy,qz),
	                                   P(rx,ry,rz), P(tx,ty,tz));
  }
};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_COPLANAR_SIDE_OF_BOUNDED_CIRCLE_3_H
