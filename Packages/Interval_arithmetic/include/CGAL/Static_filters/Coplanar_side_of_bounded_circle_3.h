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
// file          : include/CGAL/Static_filters/Coplanar_side_of_bounded_circle_3.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_STATIC_FILTERS_COPLANAR_SIDE_OF_BOUNDED_CIRCLE_3_H
#define CGAL_STATIC_FILTERS_COPLANAR_SIDE_OF_BOUNDED_CIRCLE_3_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_exact.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Profile_counter.h>
// #include <CGAL/Static_filter_error.h> // Only used to precompute constants

CGAL_BEGIN_NAMESPACE

template <class Point>
class SF_Side_of_bounded_circle_3
{
  double _static_epsilon;

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
    F det = det4x4_by_formula(t1, t1, t1, sq,
                              t1, t1, t1, sq,
                              t1, t1, t1, sq,
                              n1, n1, n1, sq_n1); // Full det
    double err = det.error();
    std::cerr << "*** epsilon for In_circle_3 = " << err << std::endl;
    return err;
  }

  static const double epsilon;

protected:

  template < class R >
  friend class Static_filters;

  // These operations are reserved to Static_filters<>, because the context of
  // a predicate is linked to the one of the Static_filter<> it is a member of.
  SF_Side_of_bounded_circle_3(const SF_Side_of_bounded_circle_3 &) {}

  SF_Side_of_bounded_circle_3& operator=(const SF_Side_of_bounded_circle_3 &)
  {}
 
  SF_Side_of_bounded_circle_3()
  {
      _static_epsilon = HUGE_VAL;
  }

public:
  typedef Bounded_side result_type;

  void update(double dx, double dy, double dz)
  {
      double d = std::max(std::max(dx, dy), dz);
      _static_epsilon = epsilon*d*d*d*d*d*d;
  }

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
    CGAL_PROFILER(calls, "In_circle_3 calls")

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

    double det = det4x4_by_formula(ptx,pty,ptz,pt2,
                                   rtx,rty,rtz,rt2,
                                   qtx,qty,qtz,qt2,
                                   vx,vy,vz,v2);

    // Try a fully static bound first.
    if (det >  _static_epsilon) return ON_BOUNDED_SIDE;
    if (det < -_static_epsilon) return ON_UNBOUNDED_SIDE;

    CGAL_PROFILER(st_fail, "In_circle_3 static failures")

    // Compute the semi-static bound.
    double maxx = fabs(px);
    if (maxx < fabs(qx)) maxx = fabs(qx);
    if (maxx < fabs(rx)) maxx = fabs(rx);
    if (maxx < fabs(tx)) maxx = fabs(tx);
    double maxy = fabs(py);
    if (maxy < fabs(qy)) maxy = fabs(qy);
    if (maxy < fabs(ry)) maxy = fabs(ry);
    if (maxy < fabs(ty)) maxy = fabs(ty);
    double maxz = fabs(pz);
    if (maxz < fabs(qz)) maxz = fabs(qz);
    if (maxz < fabs(rz)) maxz = fabs(rz);
    if (maxz < fabs(tz)) maxz = fabs(tz);

    double d = std::max(maxx, std::max(maxy, maxz));
    double eps = epsilon*d*d*d*d*d*d;

    if (det >  eps) return ON_BOUNDED_SIDE;
    if (det < -eps) return ON_UNBOUNDED_SIDE;

    CGAL_PROFILER(fail, "In_circle_3 semi-static failures")

    typedef Simple_cartesian<Filtered_exact<double, MP_Float> > K;
    typedef K::Point_3 P;

    return coplanar_side_of_bounded_circle(P(px,py,pz), P(qx,qy,qz),
	                                   P(rx,ry,rz), P(tx,ty,tz));
  }
};

template <class Point>
const double SF_Side_of_bounded_circle_3<Point>::epsilon = 3.27418e-11;

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_COPLANAR_SIDE_OF_BOUNDED_CIRCLE_3_H
