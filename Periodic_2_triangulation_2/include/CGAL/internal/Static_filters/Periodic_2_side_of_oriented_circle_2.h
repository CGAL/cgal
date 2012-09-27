// Copyright (c) 2001,2004,2008-2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://nicokruithof@scm.gforge.inria.fr/svnroot/cgal/trunk/Periodic_3_triangulation_3/include/CGAL/internal/Static_filters/Periodic_3_side_of_oriented_sphere_3.h $
// $Id: Periodic_3_side_of_oriented_sphere_3.h 63134 2011-04-26 17:01:34Z sloriot $
// 
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_2_SIDE_OF_ORIENTED_CIRCLE_2_H
#define CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_2_SIDE_OF_ORIENTED_CIRCLE_2_H

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>

#include <CGAL/Periodic_2_offset_2.h>

namespace CGAL { namespace internal { namespace Static_filters_predicates {

template < typename K_base >
class Periodic_2_side_of_oriented_circle_2
  : public K_base::Side_of_oriented_circle_2
{
  typedef typename K_base::Side_of_oriented_circle_2    Base;

  typedef typename K_base::Point_2                      Point_2;
  typedef typename K_base::Iso_rectangle_2              Iso_rectangle_2;
  typedef CGAL::Periodic_2_offset_2                     Offset;

public:
  const Iso_rectangle_2 * _dom;

public:
  typedef typename Base::result_type  result_type;

  template <class EX, class AP>
  Periodic_2_side_of_oriented_circle_2(const Iso_rectangle_2 * dom,
					    const EX * dom_e,
					    const AP * dom_f)
      : Base(dom_e,dom_f), _dom(dom) { }

  Oriented_side
  operator()(const Point_2 &p, const Point_2 &q, const Point_2 &r,
             const Point_2 &s) const
  {
      CGAL_PROFILER("Periodic_2_side_of_oriented_circle_2 calls");

      double px, py, qx, qy, rx, ry, sx, sy;

      if (fit_in_double(p.x(), px) && fit_in_double(p.y(), py) &&
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy) &&
          fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry) &&
          fit_in_double(s.x(), sx) && fit_in_double(s.y(), sy)) {

          CGAL_PROFILER("Periodic_2_side_of_oriented_circle_2 semi-static attempts");

          double psx = px - sx;
          double psy = py - sy;
          double pt2 = CGAL_NTS square(psx) + CGAL_NTS square(psy);
          double qsx = qx - sx;
          double qsy = qy - sy;
          double qt2 = CGAL_NTS square(qsx) + CGAL_NTS square(qsy);
          double rsx = rx - sx;
          double rsy = ry - sy;
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

          double eps = 1.2466136531027298e-13 * maxx * maxy * (maxy * maxy);

          double det = -CGAL::determinant(psx, psy, pt2,
                                          qsx, qsy, qt2,
                                          rsx, rsy, rt2);

          // Protect against underflow in the computation of eps.
          if (maxx < 1e-58) /* sqrt^5(min_double/eps) */ {
            if (maxx == 0)
              return ON_ORIENTED_BOUNDARY;
          }
          // Protect against overflow in the computation of det.
          else if (maxy < 1e61) /* sqrt^5(max_double/4 [hadamard]) */ {
            if (det > eps)  return ON_POSITIVE_SIDE;
            if (det < -eps) return ON_NEGATIVE_SIDE;
          }

          CGAL_PROFILER("Periodic_2_side_of_oriented_circle_2 semi-static failures");
      }
      return Base::operator()(p, q, r, s);
  }

  Oriented_side
  operator()(const Point_2 &p, const Point_2 &q,
             const Point_2 &r, const Point_2 &s,
             const Offset &o_p, const Offset &o_q,
             const Offset &o_r, const Offset &o_s) const {

    CGAL_PROFILER("Periodic_2_side_of_oriented_circle_2 calls");

    double px, py, qx, qy, rx, ry, sx, sy;
    double domxmax, domxmin, domymax, domymin;
    int osx = o_s.x();
    int osy = o_s.y();

    if (fit_in_double(p.x(), px) && fit_in_double(p.y(), py) &&
        fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy) &&
        fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry) &&
        fit_in_double(s.x(), sx) && fit_in_double(s.y(), sy) &&
	fit_in_double(_dom->xmax(), domxmax) &&
	fit_in_double(_dom->xmin(), domxmin) &&
	fit_in_double(_dom->ymax(), domymax) &&
	fit_in_double(_dom->ymin(), domymin)) {

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
      if (maxx < 1e-58) /* sqrt^5(min_double/eps) */ {
        if (maxx == 0)
          return ON_ORIENTED_BOUNDARY;
      }
      // Protect against overflow in the computation of det.
      else if (maxy < 1e61) /* sqrt^5(max_double/4 [hadamard]) */ {
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
    F t1 = F(1,F::ulp()/4);   // First translations
    F sq = t1*t1+t1*t1+t1*t1; // squares
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

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_2_SIDE_OF_ORIENTED_CIRCLE_2_H
