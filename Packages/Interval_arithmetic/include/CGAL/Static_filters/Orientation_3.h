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

#ifndef CGAL_STATIC_FILTERS_ORIENTATION_3_H
#define CGAL_STATIC_FILTERS_ORIENTATION_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Static_filter_error.h>

CGAL_BEGIN_NAMESPACE

template < typename K_base >
class SF_Orientation_3
  : public K_base::Orientation_3
{
  typedef typename K_base::Point_3          Point_3;
  typedef typename K_base::Orientation_3    Base;

public:

  Orientation operator()(const Point_3 &p, const Point_3 &q,
                         const Point_3 &r, const Point_3 &s) const
  {
      CGAL_PROFILER("Orientation_3 calls");

      double px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz;

      if (fit_in_double(p.x(), px) && fit_in_double(p.y(), py) &&
          fit_in_double(p.z(), pz) &&
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy) &&
          fit_in_double(q.z(), qz) &&
          fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry) &&
          fit_in_double(r.z(), rz) &&
          fit_in_double(s.x(), sx) && fit_in_double(s.y(), sy) &&
          fit_in_double(s.z(), sz))
      {
          CGAL_PROFILER("Orientation_3 semi-static attempts");

          double pqx = qx - px;
          double pqy = qy - py;
          double pqz = qz - pz;
          double prx = rx - px;
          double pry = ry - py;
          double prz = rz - pz;
          double psx = sx - px;
          double psy = sy - py;
          double psz = sz - pz;

          // Then semi-static filter.
          double maxx = fabs(pqx);
          if (maxx < fabs(prx)) maxx = fabs(prx);
          if (maxx < fabs(psx)) maxx = fabs(psx);
          double maxy = fabs(pqy);
          if (maxy < fabs(pry)) maxy = fabs(pry);
          if (maxy < fabs(psy)) maxy = fabs(psy);
          double maxz = fabs(pqz);
          if (maxz < fabs(prz)) maxz = fabs(prz);
          if (maxz < fabs(psz)) maxz = fabs(psz);
          double eps = 5.1107127829973276e-15 * maxx * maxy * maxz;

          double det = det3x3_by_formula(pqx, pqy, pqz,
                                         prx, pry, prz,
                                         psx, psy, psz);

          // Protect against underflow in the computation of eps.
          if (maxx < 1e-97 || maxy < 1e-97 || maxz < 1e-97) {
            if (maxx == 0 || maxy == 0 || maxz == 0)
              return ZERO;
          }
          else {
            if (det > eps)  return POSITIVE;
            if (det < -eps) return NEGATIVE;
          }

          CGAL_PROFILER("Orientation_3 semi-static failures");
      }

      return Base::operator()(p, q, r, s);
  }

  // Computes the epsilon for Orientation_3.
  static double compute_epsilon()
  {
    typedef Static_filter_error F;
    F t1 = F(1, F::ulp()/2);         // First translation
    F det = det3x3_by_formula(t1, t1, t1,
                              t1, t1, t1,
                              t1, t1, t1); // Full det
    double err = det.error();
    std::cerr << "*** epsilon for Orientation_3 = " << err << std::endl;
    return err;
  }

};

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_ORIENTATION_3_H
