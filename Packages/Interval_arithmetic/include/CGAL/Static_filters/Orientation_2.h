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

#ifndef CGAL_STATIC_FILTERS_ORIENTATION_2_H
#define CGAL_STATIC_FILTERS_ORIENTATION_2_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Profile_counter.h>
// #include <CGAL/Static_filter_error.h> // Only used to precompute constants

CGAL_BEGIN_NAMESPACE

template < typename K_base >
class SF_Orientation_2
  : public K_base::Orientation_2
{
#if 0
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
#endif

  typedef typename K_base::Point_2          Point_2;
  typedef typename K_base::Orientation_2    Base;

public:

  Orientation operator()(const Point_2 &p, const Point_2 &q, const Point_2 &r) const
  {
      if (fit_in_double(p.x()) && fit_in_double(p.y()) &&
          fit_in_double(q.x()) && fit_in_double(q.y()) &&
          fit_in_double(r.x()) && fit_in_double(r.y()))
      {
          double px = CGAL_NTS to_double(p.x());
          double py = CGAL_NTS to_double(p.y());
          double qx = CGAL_NTS to_double(q.x());
          double qy = CGAL_NTS to_double(q.y());
          double rx = CGAL_NTS to_double(r.x());
          double ry = CGAL_NTS to_double(r.y());

          CGAL_PROFILER("Orientation_2 calls");

          double pqx = qx-px;
          double pqy = qy-py;
          double prx = rx-px;
          double pry = ry-py;

          double det = det2x2_by_formula(pqx, pqy,
                                         prx, pry);

          // Then semi-static filter.
          double maxx = fabs(px);
          if (maxx < fabs(qx)) maxx = fabs(qx);
          if (maxx < fabs(rx)) maxx = fabs(rx);
          double maxy = fabs(py);
          if (maxy < fabs(qy)) maxy = fabs(qy);
          if (maxy < fabs(ry)) maxy = fabs(ry);
          double eps = 3.55271e-15 * maxx * maxy;

          if (det > eps)  return POSITIVE;
          if (det < -eps) return NEGATIVE;

          CGAL_PROFILER("Orientation_2 semi-static failures");
      }

      return Base::operator()(p, q, r);
  }

};

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_ORIENTATION_2_H
