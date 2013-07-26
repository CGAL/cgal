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
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_2_ORIENTATION_2_H
#define CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_2_ORIENTATION_2_H

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Orientation_2.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>

#include <CGAL/Periodic_2_offset_2.h>

#include <cmath>

/*
  Current implementation: Just take the 3D predicates and set z to
  zero. The epsilon bounds should still hold, although they might be
  too high.

  First attempt tried to generate the predicates using the following
  code:
  Predicates are generated using mcc from the Filtered Predicate
  Generator (fpg) using:

  double
  det2x2( double a00, double a01,
  double a10, double a11 )
  {
  return a00*a11 - a10*a01;
  }

  int
  orientationC2( double px, double py,
                 double qx, double qy,
                 double rx, double ry )
  group px qx rx;
  group py qy ry;
  {
  return sign( det2x2( qx-px, qy-py ,
  rx-px, ry-py ) );
  }

  int
  orientationC2_offset( double px, double py, int opx, int opy,
                        double qx, double qy, int oqx, int oqy,
                        double rx, double ry, int orx, int ory,
                        double domxmax, double domxmin,
                        double domymax, double domymin
  )
  group px qx rx domxmax domxmin opx oqx orx;
  group py qy ry domymax domymin opy oqy ory;
  {
  double domx = domxmax-domxmin;
  double domy = domymax-domymin;

  double pqx = 1*(qx - px) + (oqx-opx)*domx;
  double pqy = 1*(qy - py) + (oqy-opy)*domy;
  double prx = 1*(rx - px) + (orx-opx)*domx;
  double pry = 1*(ry - py) + (ory-opy)*domy;

  return sign( det2x2(pqx, pqy, prx, pry) );
  }
*/

namespace CGAL
{
namespace internal
{
namespace Static_filters_predicates
{

template < typename K_base >
class Periodic_2_orientation_2 : public K_base::Orientation_2
{
  typedef typename K_base::Point_2          Point_2;
  typedef typename K_base::Vector_2         Vector_2;
  typedef typename K_base::Circle_2         Circle_2;
  typedef typename K_base::Offset_2         Offset_2;
  typedef typename K_base::Iso_rectangle_2  Iso_rectangle_2;

  typedef typename K_base::Orientation_2    Base;

public:
  const Iso_rectangle_2 * const _dom;

public:
  typedef typename Base::result_type  result_type;

  template <class EX, class AP>
  Periodic_2_orientation_2(const Iso_rectangle_2 * const dom,
                           const EX * dom_e, const AP * dom_f) : Base(dom_e, dom_f), _dom(dom)
  {
  }

#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else
  result_type
  operator()(const Vector_2& u, const Vector_2& v) const
  {
    return Base::operator()(u, v);
  }

  result_type
  operator()(const Circle_2& c) const
  {
    return Base::operator()(c);
  }
#endif

  /// Normal static orientation test, copied from Orientation_2
  result_type operator()(const Point_2 &p, const Point_2 &q, const Point_2 &r) const
  {
    CGAL_PROFILER("Periodic_2_orientation_2 calls");

    double px, py, qx, qy, rx, ry;

    if (fit_in_double(p.x(), px) && fit_in_double(p.y(), py) &&
        fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy) &&
        fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry))
      {
        CGAL_PROFILER("Periodic_2_orientation_2 semi-static attempts");

        double pqx = qx - px;
        double pqy = qy - py;
        double prx = rx - px;
        double pry = ry - py;

        // Then semi-static filter.
        double maxx = CGAL::abs(pqx);
        double maxy = CGAL::abs(pqy);

        double aprx = CGAL::abs(prx);
        double apry = CGAL::abs(pry);

        if (maxx < aprx) maxx = aprx;
        if (maxy < apry) maxy = apry;
        double eps = 5.1107127829973299e-15 * maxx * maxy;
        double det = CGAL::determinant(pqx, pqy,
                                       prx, pry);

        // Sort maxx < maxy
        if (maxx > maxy)
          std::swap(maxx, maxy);

        // Protect against underflow in the computation of eps.
        if (maxx < 1e-97) /* cbrt(min_double/eps) */
          {
            if (maxx == 0)
              return ZERO;
          }
        // Protect against overflow in the computation of det.
        else if (maxy < 1e102) /* cbrt(max_double [hadamard]/4) */
          {
            if (det > eps)  return POSITIVE;
            if (det < -eps) return NEGATIVE;
          }

        CGAL_PROFILER("Periodic_2_orientation_2 semi-static failures");
      }

    return Base::operator()(p, q, r);
  }


  /// Static orientation test with offsets
  result_type operator()(const Point_2 &p, const Point_2 &q, const Point_2 &r,
                         const Offset_2 &o_p, const Offset_2 &o_q, const Offset_2 &o_r) const
  {

    CGAL_PROFILER("Periodic_2_orientation_2 with offset calls");

    double px, py, qx, qy, rx, ry;
    double domxmax, domxmin, domymax, domymin;
    int opx = o_p.x();
    int opy = o_p.y();

    if (fit_in_double(p.x(), px) && fit_in_double(p.y(), py) &&
        fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy) &&
        fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry) &&
        fit_in_double(_dom->xmax(), domxmax) &&
        fit_in_double(_dom->xmin(), domxmin) &&
        fit_in_double(_dom->ymax(), domymax) &&
        fit_in_double(_dom->ymin(), domymin))
      {
        CGAL_PROFILER("Periodic_2_orientation_2 with offset semi-static attempts");

        double domx = domxmax - domxmin;
        double domy = domymax - domymin;

        double pqx = qx - px + domx * ( o_q.x() - opx );
        double pqy = qy - py + domy * ( o_q.y() - opy );
        double prx = rx - px + domx * ( o_r.x() - opx );
        double pry = ry - py + domy * ( o_r.y() - opy );

        // Then semi-static filter.
        double maxx = CGAL::abs(pqx);
        double maxy = CGAL::abs(pqy);

        double aprx = CGAL::abs(prx);
        double apry = CGAL::abs(pry);

        if (maxx < aprx) maxx = aprx;
        if (maxy < apry) maxy = apry;
        double eps = 4.111024169857068197e-15 * maxx * maxy;
        double det = CGAL::determinant(pqx, pqy,
                                       prx, pry);

        // Sort maxx < maxy.
        if (maxx > maxy)
          std::swap(maxx, maxy);

        // Protect against underflow in the computation of eps.
        if (maxx < 1e-97) /* cbrt(min_double/eps) */
          {
            if (maxx == 0)
              return ZERO;
          }
        // Protect against overflow in the computation of det.
        else if (maxy < 1e102) /* cbrt(max_double [hadamard]/4) */
          {
            if (det > eps)  return POSITIVE;
            if (det < -eps) return NEGATIVE;
          }

        CGAL_PROFILER("Periodic_2_orientation_2 with offset semi-static failures");
      }

    return Base::operator()(p, q, r, o_p, o_q, o_r);
  }

  // Computes the epsilon for Periodic_2_orientation_2.
  static double compute_epsilon()
  {
    typedef Static_filter_error F;
    F t1 = F(1, F::ulp() / 4);       // First translation
    F det = CGAL::determinant(t1, t1, t1,
                              t1, t1, t1,
                              t1, t1, t1); // Full det
    double err = det.error();
    err += err * 2 * F::ulp(); // Correction due to "eps * maxx * maxy...".
    std::cerr << "*** epsilon for Periodic_2_orientation_2 = " << err
              << std::endl;
    return err;
  }
};

}
}
} // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_2_ORIENTATION_2_H
