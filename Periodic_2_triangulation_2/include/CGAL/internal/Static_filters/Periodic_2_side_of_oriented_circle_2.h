// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_2_SIDE_OF_ORIENTED_CIRCLE_2_H
#define CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_2_SIDE_OF_ORIENTED_CIRCLE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/internal/Static_filters/tools.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <CGAL/internal/Static_filters/tools.h>

#include <CGAL/Periodic_2_offset_2.h>

namespace CGAL {
namespace internal {
namespace Static_filters_predicates {

template <class K_, class Side_of_oriented_circle_2_base_>
class Periodic_2_side_of_oriented_circle_2
  : public Side_of_oriented_circle_2_base_
{
  typedef Side_of_oriented_circle_2_base_               Base;

public:
  typedef K_                                            Kernel;

  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Iso_rectangle_2              Iso_rectangle_2;
  typedef typename Kernel::Periodic_2_offset_2          Offset;

private:
  const Iso_rectangle_2 * _dom;

public:
  typedef typename Base::result_type  result_type;

  Periodic_2_side_of_oriented_circle_2(const Iso_rectangle_2 * dom,
                                       const Base& socb)
    : Base(socb), _dom(dom)
  { }

  Oriented_side
  operator()(const Point_2 &p, const Point_2 &q, const Point_2 &r,
             const Point_2 &t) const
  {
    CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Periodic_2_side_of_oriented_circle_2", tmp);

    Get_approx<Point_2> get_approx; // Identity functor for all points
    // but lazy points.

    double px, py, qx, qy, rx, ry, tx, ty;
    init_double(px, py, qx, qy, rx, ry, tx, ty, (FT*)(0));
    if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
        fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
        fit_in_double(get_approx(r).x(), rx) && fit_in_double(get_approx(r).y(), ry) &&
        fit_in_double(get_approx(t).x(), tx) && fit_in_double(get_approx(t).y(), ty))
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);

        double qpx = qx - px;
        double qpy = qy - py;
        double rpx = rx - px;
        double rpy = ry - py;
        double tpx = tx - px;
        double tpy = ty - py;

        double tqx = tx - qx;
        double tqy = ty - qy;
        double rqx = rx - qx;
        double rqy = ry - qy;

        double det = CGAL::determinant(qpx * tpy - qpy * tpx, tpx * tqx + tpy * tqy,
                                       qpx * rpy - qpy * rpx, rpx * rqx + rpy * rqy);

        // We compute the semi-static bound.
        double maxx = CGAL::abs(qpx);
        double maxy = CGAL::abs(qpy);

        double arpx = CGAL::abs(rpx);
        double arpy = CGAL::abs(rpy);

        double atqx = CGAL::abs(tqx);
        double atqy = CGAL::abs(tqy);

        double atpx = CGAL::abs(tpx);
        double atpy = CGAL::abs(tpy);

        double arqx = CGAL::abs(rqx);
        double arqy = CGAL::abs(rqy);

        if (maxx < arpx) maxx = arpx;
        if (maxx < atpx) maxx = atpx;
        if (maxx < atqx) maxx = atqx;
        if (maxx < arqx) maxx = arqx;

        if (maxy < arpy) maxy = arpy;
        if (maxy < atpy) maxy = atpy;
        if (maxy < atqy) maxy = atqy;
        if (maxy < arqy) maxy = arqy;

        if (maxx > maxy)  std::swap(maxx, maxy);

        // Protect against underflow in the computation of eps.
        if (maxx < 1e-73)
          {
            if (maxx == 0)
              return ON_ORIENTED_BOUNDARY;
          }
        else if (maxy < 1e76) /* sqrt(sqrt(max_double/16 [hadamard])) */
          {
            double eps = 8.8878565762001373e-15 * maxx * maxy * (maxy * maxy);
            if (det > eps)  return ON_POSITIVE_SIDE;
            if (det < -eps) return ON_NEGATIVE_SIDE;
          }

        CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
      }

    return Base::operator()(p, q, r, t);
  }

  Oriented_side
  operator()(const Point_2 &p, const Point_2 &q,
             const Point_2 &r, const Point_2 &s,
             const Offset &o_p, const Offset &o_q,
             const Offset &o_r, const Offset &o_s) const
  {

    CGAL_PROFILER("Periodic_2_side_of_oriented_circle_2 calls");

    double px, py, qx, qy, rx, ry, sx, sy;
    double domxmax, domxmin, domymax, domymin;
    int osx = o_s.x();
    int osy = o_s.y();
    init_double(px, py, qx, qy, rx, ry, sx, sy, domxmax, domxmin, domymax, domymin, (FT*)(0));
    if (fit_in_double(p.x(), px) && fit_in_double(p.y(), py) &&
        fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy) &&
        fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry) &&
        fit_in_double(s.x(), sx) && fit_in_double(s.y(), sy) &&
        fit_in_double(_dom->xmax(), domxmax) &&
        fit_in_double(_dom->xmin(), domxmin) &&
        fit_in_double(_dom->ymax(), domymax) &&
        fit_in_double(_dom->ymin(), domymin))
      {

        CGAL_PROFILER("Periodic_2_side_of_oriented_circle_2 with offset semi-static attempts");

        double domx = domxmax - domxmin;
        double domy = domymax - domymin;

        double psx = px - sx + domx * (o_p.x() - osx);
        double psy = py - sy + domy * (o_p.y() - osy);
        double pt2 = CGAL_NTS square(psx) + CGAL_NTS square(psy);
        double qsx = qx - sx + domx * (o_q.x() - osx);
        double qsy = qy - sy + domy * (o_q.y() - osy);
        double qt2 = CGAL_NTS square(qsx) + CGAL_NTS square(qsy);
        double rsx = rx - sx + domx * (o_r.x() - osx);
        double rsy = ry - sy + domy * (o_r.y() - osy);
        double rt2 = CGAL_NTS square(rsx) + CGAL_NTS square(rsy);

        // Compute the semi-static bound.
        double maxx = CGAL::abs(psx);
        double maxy = CGAL::abs(psy);

        double aqsx = CGAL::abs(qsx);
        double aqsy = CGAL::abs(qsy);

        double arsx = CGAL::abs(rsx);
        double arsy = CGAL::abs(rsy);

        if (maxx < aqsx) maxx = aqsx;
        if (maxx < arsx) maxx = arsx;

        if (maxy < aqsy) maxy = aqsy;
        if (maxy < arsy) maxy = arsy;

        // Sort maxx < maxy.
        if (maxx > maxy)
          std::swap(maxx, maxy);

        double eps = 1.0466759304746772485e-13 * maxx * maxy * (maxy * maxy);

        double det = CGAL::determinant(psx, psy, pt2,
                                       qsx, qsy, qt2,
                                       rsx, rsy, rt2);

        // Protect against underflow in the computation of eps.
        if (maxx < 1e-58) /* sqrt^5(min_double/eps) */
          {
            if (maxx == 0)
              return ON_ORIENTED_BOUNDARY;
          }
        // Protect against overflow in the computation of det.
        else if (maxy < 1e61) /* sqrt^5(max_double/4 [hadamard]) */
          {
            if (det > eps)  return ON_POSITIVE_SIDE;
            if (det < -eps) return ON_NEGATIVE_SIDE;
          }

        CGAL_PROFILER("Periodic_2_side_of_oriented_circle_2 with offset semi-static failures");
      }
    return Base::operator()(p, q, r, s, o_p, o_q, o_r, o_s);
  }

  // Computes the epsilon for Periodic_2_side_of_oriented_circle_2.
  static double compute_epsilon()
  {
    typedef Static_filter_error F;
    F t1 = F(1, F::ulp() / 4); // First translations
    F sq = t1 * t1 + t1 * t1 + t1 * t1; // squares
    F det = CGAL::determinant(t1, t1, t1, sq,
                              t1, t1, t1, sq,
                              t1, t1, t1, sq,
                              t1, t1, t1, sq); // Full det
    double err = det.error();
    err += err * 3 * F::ulp(); // Correction due to "eps * maxx * ...".

    std::cerr << "*** epsilon for Periodic_2_side_of_oriented_circle_2 = "
              << err << std::endl;
    return err;
  }
};

}
}
} // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_2_SIDE_OF_ORIENTED_CIRCLE_2_H
