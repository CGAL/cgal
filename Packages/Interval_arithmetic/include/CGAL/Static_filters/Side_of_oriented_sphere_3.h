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

#ifndef CGAL_STATIC_FILTERS_SIDE_OF_ORIENTED_SPHERE_3_H
#define CGAL_STATIC_FILTERS_SIDE_OF_ORIENTED_SPHERE_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Static_filter_error.h>

CGAL_BEGIN_NAMESPACE

template < typename K_base >
class SF_Side_of_oriented_sphere_3
  : public K_base::Side_of_oriented_sphere_3
{
  typedef typename K_base::Point_3                      Point_3;
  typedef typename K_base::Side_of_oriented_sphere_3    Base;

public:

  Oriented_side
  operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r,
             const Point_3 &s, const Point_3 &t) const
  {
      CGAL_PROFILER("In_sphere_3 calls");

      if (fit_in_double(p.x()) && fit_in_double(p.y()) && fit_in_double(p.z()) &&
          fit_in_double(q.x()) && fit_in_double(q.y()) && fit_in_double(q.z()) &&
          fit_in_double(r.x()) && fit_in_double(r.y()) && fit_in_double(r.z()) &&
          fit_in_double(s.x()) && fit_in_double(s.y()) && fit_in_double(s.z()) &&
          fit_in_double(t.x()) && fit_in_double(t.y()) && fit_in_double(t.z()))
      {
          CGAL_PROFILER("In_sphere_3 semi-static attempts");

          const double & tx = CGAL_NTS to_double(t.x());
          const double & ty = CGAL_NTS to_double(t.y());
          const double & tz = CGAL_NTS to_double(t.z());

          double ptx = CGAL_NTS to_double(p.x()) - tx;
          double pty = CGAL_NTS to_double(p.y()) - ty;
          double ptz = CGAL_NTS to_double(p.z()) - tz;
          double pt2 = CGAL_NTS square(ptx) + CGAL_NTS square(pty)
	             + CGAL_NTS square(ptz);
          double qtx = CGAL_NTS to_double(q.x()) - tx;
          double qty = CGAL_NTS to_double(q.y()) - ty;
          double qtz = CGAL_NTS to_double(q.z()) - tz;
          double qt2 = CGAL_NTS square(qtx) + CGAL_NTS square(qty)
	             + CGAL_NTS square(qtz);
          double rtx = CGAL_NTS to_double(r.x()) - tx;
          double rty = CGAL_NTS to_double(r.y()) - ty;
          double rtz = CGAL_NTS to_double(r.z()) - tz;
          double rt2 = CGAL_NTS square(rtx) + CGAL_NTS square(rty)
	             + CGAL_NTS square(rtz);
          double stx = CGAL_NTS to_double(s.x()) - tx;
          double sty = CGAL_NTS to_double(s.y()) - ty;
          double stz = CGAL_NTS to_double(s.z()) - tz;
          double st2 = CGAL_NTS square(stx) + CGAL_NTS square(sty)
	             + CGAL_NTS square(stz);

          // Compute the semi-static bound.
          double maxx = fabs(ptx);
          if (maxx < fabs(qtx)) maxx = fabs(qtx);
          if (maxx < fabs(rtx)) maxx = fabs(rtx);
          if (maxx < fabs(stx)) maxx = fabs(stx);
          double maxy = fabs(pty);
          if (maxy < fabs(qty)) maxy = fabs(qty);
          if (maxy < fabs(rty)) maxy = fabs(rty);
          if (maxy < fabs(sty)) maxy = fabs(sty);
          double maxz = fabs(ptz);
          if (maxz < fabs(qtz)) maxz = fabs(qtz);
          if (maxz < fabs(rtz)) maxz = fabs(rtz);
          if (maxz < fabs(stz)) maxz = fabs(stz);

	  double maxt = maxx;
          if (maxt < maxy) maxt = maxy;
          if (maxt < maxz) maxt = maxz;

          double eps = 1.14631e-13 * maxx * maxy * maxz * (maxt * maxt);

          double det = det4x4_by_formula(ptx,pty,ptz,pt2,
                                         rtx,rty,rtz,rt2,
                                         qtx,qty,qtz,qt2,
                                         stx,sty,stz,st2);

          if (det >  eps) return ON_POSITIVE_SIDE;
          if (det < -eps) return ON_NEGATIVE_SIDE;

          CGAL_PROFILER("In_sphere_3 semi-static failures");
      }
      return Base::operator()(p, q, r, s, t);
  }

  // Computes the epsilon for In_sphere_3.
  static double compute_epsilon()
  {
    typedef CGAL::Static_filter_error F;
    F t1 = F(1,F::ulp()/2);         // First translation
    F sq = t1*t1+t1*t1+t1*t1; // squares
    F det = det4x4_by_formula(t1, t1, t1, sq,
                              t1, t1, t1, sq,
                              t1, t1, t1, sq,
                              t1, t1, t1, sq); // Full det
    double err = det.error();
    std::cerr << "*** epsilon for In_sphere_3 = " << err << std::endl;
    return err;
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_SIDE_OF_ORIENTED_SPHERE_3_H
