// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Static_filters/Side_of_oriented_circle_2.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

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
  double _static_epsilon;

public:
  // Computes the epsilon for In_circle_2.
  static void cir_2()
  {
    Static_filter_error X(1);
    side_of_oriented_circleC2(X, X, X, X, X, X, X, X);
  }

  static const double epsilon; // = 1.42109e-13; // cir_2();

protected:

  template < class R >
  friend class Static_filters;

  // These operations are reserved to Static_filters<>, because the context of
  // a predicate is linked to the one of the Static_filter<> it is a member of.
  SF_Side_of_oriented_circle_2(const SF_Side_of_oriented_circle_2 &s)
      : _static_epsilon(s._static_epsilon) {}

  SF_Side_of_oriented_circle_2& operator=(const SF_Side_of_oriented_circle_2&s)
  {
      _static_epsilon = s._static_epsilon;
      return *this;
  }
 
  SF_Side_of_oriented_circle_2()
  {
      _static_epsilon = HUGE_VAL;
  }

public:
  typedef Oriented_side result_type;

  void update(double dx, double dy)
  {
      _static_epsilon = epsilon*dx*dy*2*(dx*dx+dy*dy);
  }

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
    CGAL_PROFILER(calls, "In_circle_2 calls")

    double qpx = qx-px;
    double qpy = qy-py;
    double rpx = rx-px;
    double rpy = ry-py;
    double tpx = tx-px;
    double tpy = ty-py;

    double det = det2x2_by_formula(
                             qpx*tpy - qpy*tpx, tpx*(tx-qx) + tpy*(ty-qy),
                             qpx*rpy - qpy*rpx, rpx*(rx-qx) + rpy*(ry-qy));

    // Try a fully static bound first, when possible.
    if (det >  _static_epsilon) return ON_POSITIVE_SIDE;
    if (det < -_static_epsilon) return ON_NEGATIVE_SIDE;

    CGAL_PROFILER(st_fail, "In_circle_2 static failures")

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

    CGAL_PROFILER(fail, "In_circle_2 semi-static failures")

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
	CGAL_PROFILER(exact_diff, "In_circle_2 exact diffs")

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

	CGAL_PROFILER(step2, "In_circle_2 step2 failures")
    }

    CGAL_PROFILER(step3, "In_circle_2 step3")

    typedef Simple_cartesian<Filtered_exact<double, MP_Float> > K;
    typedef K::Point_2 P;

    Oriented_side oooo = side_of_oriented_circle(P(px,py), P(qx,qy),
	                                         P(rx,ry), P(tx,ty));
    if (oooo == ON_ORIENTED_BOUNDARY) {
	CGAL_PROFILER(is_null, "In_circle_2 is_null")
    }
    return oooo;
  }
};

template <class Point>
const double SF_Side_of_oriented_circle_2<Point>::epsilon = 1.42109e-13;

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_SIDE_OF_ORIENTED_CIRCLE_2_H
