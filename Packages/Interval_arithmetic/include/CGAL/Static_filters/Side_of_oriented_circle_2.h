// Copyright (c) 2001,2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_STATIC_FILTERS_SIDE_OF_ORIENTED_CIRCLE_2_H
#define CGAL_STATIC_FILTERS_SIDE_OF_ORIENTED_CIRCLE_2_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_exact.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Profile_counter.h>
// #include <CGAL/Static_filter_error.h> // Only used to precompute constants

CGAL_BEGIN_NAMESPACE

template <class Point>
class SF_Side_of_oriented_circle_2
{
public:
  // Computes the epsilon for In_circle_2.
  static void cir_2()
  {
    Static_filter_error X(1);
    side_of_oriented_circleC2(X, X, X, X, X, X, X, X);
  }

  static const double epsilon; // = 1.42109e-13; // cir_2();

public:
  typedef Oriented_side result_type;

  Oriented_side operator()(const Point &p, const Point &q,
	                   const Point &r, const Point &t) const
  {
    return opti_in_circleC2(
	    to_double(p.x()), to_double(p.y()),
            to_double(q.x()), to_double(q.y()),
	    to_double(r.x()), to_double(r.y()),
	    to_double(t.x()), to_double(t.y()));
  }

  Oriented_side
  opti_in_circleC2(double px, double py,
                   double qx, double qy,
		   double rx, double ry,
		   double tx, double ty) const
  {
    CGAL_PROFILER("In_circle_2 calls");

    double qpx = qx-px;
    double qpy = qy-py;
    double rpx = rx-px;
    double rpy = ry-py;
    double tpx = tx-px;
    double tpy = ty-py;

    double det = det2x2_by_formula(
                             qpx*tpy - qpy*tpx, tpx*(tx-qx) + tpy*(ty-qy),
                             qpx*rpy - qpy*rpx, rpx*(rx-qx) + rpy*(ry-qy));

    // We compute the semi-static bound.
    double maxx = fabs(px);
    if (maxx < fabs(qx)) maxx = fabs(qx);
    if (maxx < fabs(rx)) maxx = fabs(rx);
    if (maxx < fabs(tx)) maxx = fabs(tx);
    double maxy = fabs(py);
    if (maxy < fabs(qy)) maxy = fabs(qy);
    if (maxy < fabs(ry)) maxy = fabs(ry);
    if (maxy < fabs(ty)) maxy = fabs(ty);

    double pp = px*px + py*py;
    double qq = qx*qx + qy*qy;
    double rr = rx*rx + ry*ry;
    double max2 = tx*tx + ty*ty;
    if (max2 < qq) max2 = qq;
    if (max2 < rr) max2 = rr;
    double eps = epsilon*maxx*maxy*(max2+pp);

    if (det >  eps) return ON_POSITIVE_SIDE;
    if (det < -eps) return ON_NEGATIVE_SIDE;

    CGAL_PROFILER("In_circle_2 semi-static failures");

    // This predicate is different from Orientation in that all arguments are
    // local.  Thus the differences have a big probability to have been exact,
    // and helps a lot in reducing the bound of the last column.

    if (diff_was_exact(qx, px, qpx) &&
        diff_was_exact(qy, py, qpy) &&
        diff_was_exact(rx, px, rpx) &&
        diff_was_exact(ry, py, rpy) &&
        diff_was_exact(tx, px, tpx) &&
        diff_was_exact(ty, py, tpy))
    {
	CGAL_PROFILER("In_circle_2 exact diffs");

        double max2 = tpx*tpx + tpy*tpy;
        double qq = qpx*qpx + qpy*qpy;
        double rr = rpx*rpx + rpy*rpy;
        if (max2 < qq) max2 = qq;
        if (max2 < rr) max2 = rr;
        // maxx, maxy can be based on ptx and co directly, now...
        double maxx = fabs(tpx);
        if (maxx < fabs(qpx)) maxx = fabs(qpx);
        if (maxx < fabs(rpx)) maxx = fabs(rpx);
        double maxy = fabs(tpy);
        if (maxy < fabs(qpy)) maxy = fabs(qpy);
        if (maxy < fabs(rpy)) maxy = fabs(rpy);
        double eps = epsilon*maxx*maxy*max2;

        if (det >  eps) return ON_POSITIVE_SIDE;
        if (det < -eps) return ON_NEGATIVE_SIDE;

	CGAL_PROFILER("In_circle_2 step2 failures");
    }

    CGAL_PROFILER("In_circle_2 step3");

    typedef Simple_cartesian<Filtered_exact<double, MP_Float> > K;
    typedef K::Point_2 P;

    Oriented_side oooo = side_of_oriented_circle(P(px,py), P(qx,qy),
	                                         P(rx,ry), P(tx,ty));
    if (oooo == ON_ORIENTED_BOUNDARY) {
	CGAL_PROFILER("In_circle_2 is_null");
    }
    return oooo;
  }
};

template <class Point>
const double SF_Side_of_oriented_circle_2<Point>::epsilon = 1.42109e-13;

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_SIDE_OF_ORIENTED_CIRCLE_2_H
