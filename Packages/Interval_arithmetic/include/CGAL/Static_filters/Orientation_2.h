// Copyright (c) 2001  Utrecht University (The Netherlands),
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

#ifndef CGAL_STATIC_FILTERS_ORIENTATION_2_H
#define CGAL_STATIC_FILTERS_ORIENTATION_2_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_exact.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Profile_counter.h>
// #include <CGAL/Static_filter_error.h> // Only used to precompute constants

CGAL_BEGIN_NAMESPACE

template < class R > class Static_filters;
template < class Pt, class Or2 > class SF_Coplanar_orientation_3;

template <class Point>
class SF_Orientation_2
{
  double _static_epsilon;

  // Computes the epsilon for Orientation_2.
  static double ori_2()
  {
    typedef Static_filter_error F;
    F t1 = F(1)-F(1);         // First translation
    F det = det2x2_by_formula(t1, t1,
                              t1, t1); // Full det
    double err = det.error();
    std::cerr << "*** epsilon for Orientation_2 = " << err << std::endl;
    return err;
  }

  static const double epsilon; // = 3.55271e-15; // ori_2();

protected:

  template < class R >
  friend class Static_filters;

  template < class Pt, class Or2 >
  friend class SF_Coplanar_orientation_3;

  // These operations are reserved to Static_filters<>, because the context of
  // a predicate is linked to the one of the Static_filter<> it is a member of.
  SF_Orientation_2(const SF_Orientation_2 &s)
      : _static_epsilon(s._static_epsilon) {}

  SF_Orientation_2 & operator=(const SF_Orientation_2 &s)
  {
      _static_epsilon = s._static_epsilon;
      return *this;
  }

  SF_Orientation_2()
  {
      _static_epsilon = CGALi::infinity;
  }

public:
  typedef Orientation result_type;

  void update(double dx, double dy)
  {
      _static_epsilon = dx*dy*epsilon;
  }

  Orientation operator()(const Point &p, const Point &q, const Point &r) const
  {
    return opti_orientationC2(to_double(p.x()), to_double(p.y()),
	                      to_double(q.x()), to_double(q.y()),
	                      to_double(r.x()), to_double(r.y()));
  }

private:

  typedef Simple_cartesian<Filtered_exact<double, MP_Float> >::Point_2 P;

  Orientation
  opti_orientationC2(double px, double py,
                     double qx, double qy,
		     double rx, double ry) const
  {
    CGAL_PROFILER("Orientation_2 calls");

    double pqx = qx-px;
    double pqy = qy-py;
    double prx = rx-px;
    double pry = ry-py;

    double det = det2x2_by_formula(pqx, pqy,
                                   prx, pry);

    // Fully static filter first.
    if (det >  _static_epsilon) return POSITIVE;
    if (det < -_static_epsilon) return NEGATIVE;

    CGAL_PROFILER("Orientation_2 static failures");

    // Then semi-static filter.
    double maxx = fabs(px);
    if (maxx < fabs(qx)) maxx = fabs(qx);
    if (maxx < fabs(rx)) maxx = fabs(rx);
    double maxy = fabs(py);
    if (maxy < fabs(qy)) maxy = fabs(qy);
    if (maxy < fabs(ry)) maxy = fabs(ry);
    double eps = epsilon*maxx*maxy;

    if (det > eps)  return POSITIVE;
    if (det < -eps) return NEGATIVE;

    CGAL_PROFILER("Orientation_2 semi-static failures");

    Orientation oooo = orientation(P(px,py), P(qx,qy), P(rx,ry));
    if (oooo == ZERO) {
	CGAL_PROFILER("Orientation_2 det_is_null");
    }
    return oooo;
  }

};

template <class Point>
const double SF_Orientation_2<Point>::epsilon = 3.55271e-15;

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_ORIENTATION_2_H
