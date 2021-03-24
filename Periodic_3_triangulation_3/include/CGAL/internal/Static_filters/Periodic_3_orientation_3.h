// Copyright (c) 2001,2004,2008-2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_3_ORIENTATION_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_3_ORIENTATION_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <CGAL/internal/Static_filters/tools.h>

#include <CGAL/Periodic_3_offset_3.h>

#include <cmath>

namespace CGAL { namespace internal { namespace Static_filters_predicates {

template < class K, class Orientation_3_base >
class Periodic_3_orientation_3
  : public Orientation_3_base
{
  typedef Orientation_3_base                 Base;

public:
  typedef K                                  Kernel;

  typedef typename K::FT                     FT;
  typedef typename K::Point_3                Point_3;
  typedef typename K::Vector_3               Vector_3;
  typedef typename K::Iso_cuboid_3           Iso_cuboid_3;
  typedef typename K::Sphere_3               Sphere_3;
  typedef typename K::Periodic_3_offset_3    Offset;

public:
  const Iso_cuboid_3 * const _dom;

public:
 typedef typename Base::result_type  result_type;

 Periodic_3_orientation_3(const Iso_cuboid_3 * const dom,
                          const Orientation_3_base& o3b)
   : Base(o3b), _dom(dom)
 { }

  using Base::operator();

  result_type
  operator()(const Point_3 &p, const Point_3 &q,
             const Point_3 &r, const Point_3 &s) const
  {
      CGAL_PROFILER("Periodic_3_orientation_3 calls");
      Get_approx<Point_3> get_approx; // Identity functor for all points
                                      // but lazy points.
      double px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz;
      init_double(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, (FT*)(0));

      if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
          fit_in_double(get_approx(p).z(), pz) &&
          fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
          fit_in_double(get_approx(q).z(), qz) &&
          fit_in_double(get_approx(r).x(), rx) && fit_in_double(get_approx(r).y(), ry) &&
          fit_in_double(get_approx(r).z(), rz) &&
          fit_in_double(get_approx(s).x(), sx) && fit_in_double(get_approx(s).y(), sy) &&
          fit_in_double(get_approx(s).z(), sz))
      {
          CGAL_PROFILER("Periodic_3_orientation_3 semi-static attempts");

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
          double maxx = CGAL::abs(pqx);
          double maxy = CGAL::abs(pqy);
          double maxz = CGAL::abs(pqz);

          double aprx = CGAL::abs(prx);
          double apsx = CGAL::abs(psx);

          double apry = CGAL::abs(pry);
          double apsy = CGAL::abs(psy);

          double aprz = CGAL::abs(prz);
          double apsz = CGAL::abs(psz);

#ifdef CGAL_USE_SSE2_MAX
          CGAL::Max<double> mmax;

          maxx = mmax(maxx, aprx, apsx);
          maxy = mmax(maxy, apry, apsy);
          maxz = mmax(maxz, aprz, apsz);
#else
          if (maxx < aprx) maxx = aprx;
          if (maxx < apsx) maxx = apsx;
          if (maxy < apry) maxy = apry;
          if (maxy < apsy) maxy = apsy;
          if (maxz < aprz) maxz = aprz;
          if (maxz < apsz) maxz = apsz;
#endif

          double det = CGAL::determinant(pqx, pqy, pqz,
                                         prx, pry, prz,
                                         psx, psy, psz);

          double eps = 5.1107127829973299e-15 * maxx * maxy * maxz;

#ifdef CGAL_USE_SSE2_MAX
          /*
          CGAL::Min<double> mmin;
          double tmp = mmin(maxx, maxy, maxz);
          maxz = mmax(maxx, maxy, maxz);
          maxx = tmp;
          */
          sse2minmax(maxx,maxy,maxz);
          // maxy can contain ANY element
#else
          // Sort maxx < maxy < maxz.
          if (maxx > maxz)
              std::swap(maxx, maxz);
          if (maxy > maxz)
              std::swap(maxy, maxz);
          else if (maxy < maxx)
              std::swap(maxx, maxy);
#endif

          // Protect against underflow in the computation of eps.
          if (maxx < 1e-97) /* cbrt(min_double/eps) */ {
            if (maxx == 0)
              return ZERO;
          }
          // Protect against overflow in the computation of det.
          else if (maxz < 1e102) /* cbrt(max_double [hadamard]/4) */ {
            if (det > eps)  return POSITIVE;
            if (det < -eps) return NEGATIVE;
          }

          CGAL_PROFILER("Periodic_3_orientation_3 semi-static failures");
      }

      return Base::operator()(p, q, r, s);
  }

  result_type
  operator()(const Point_3 &p, const Point_3 &q,
             const Point_3 &r, const Point_3 &s,
             const Offset &o_p, const Offset &o_q,
             const Offset &o_r, const Offset &o_s) const
  {
      CGAL_PROFILER("Periodic_3_orientation_3 calls");
      Get_approx<Point_3> get_approx; // Identity functor for all points
                                      // but lazy points.
      double px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz;
      double domxmax, domxmin, domymax, domymin, domzmax, domzmin;
      init_double(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, (FT*)(0));
      init_double(domxmax, domxmin, domymax, domymin, domzmax, domzmin, (FT*)(0));
      int opx = o_p.x();
      int opy = o_p.y();
      int opz = o_p.z();

      if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
          fit_in_double(get_approx(p).z(), pz) &&
          fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
          fit_in_double(get_approx(q).z(), qz) &&
          fit_in_double(get_approx(r).x(), rx) && fit_in_double(get_approx(r).y(), ry) &&
          fit_in_double(get_approx(r).z(), rz) &&
          fit_in_double(get_approx(s).x(), sx) && fit_in_double(get_approx(s).y(), sy) &&
          fit_in_double(get_approx(s).z(), sz) &&
          fit_in_double(_dom->xmax(), domxmax) &&
          fit_in_double(_dom->xmin(), domxmin) &&
          fit_in_double(_dom->ymax(), domymax) &&
          fit_in_double(_dom->ymin(), domymin) &&
          fit_in_double(_dom->zmax(), domzmax) &&
          fit_in_double(_dom->zmin(), domzmin))
      {
          CGAL_PROFILER("Periodic_3_orientation_3 semi-static attempts");

          double domx = domxmax - domxmin;
          double domy = domymax - domymin;
          double domz = domzmax - domzmin;

          double pqx = qx - px + domx * ( o_q.x() - opx );
          double pqy = qy - py + domy * ( o_q.y() - opy );
          double pqz = qz - pz + domz * ( o_q.z() - opz );
          double prx = rx - px + domx * ( o_r.x() - opx );
          double pry = ry - py + domy * ( o_r.y() - opy );
          double prz = rz - pz + domz * ( o_r.z() - opz );
          double psx = sx - px + domx * ( o_s.x() - opx );
          double psy = sy - py + domy * ( o_s.y() - opy );
          double psz = sz - pz + domz * ( o_s.z() - opz );

          // Then semi-static filter.
          double maxx = CGAL::abs(pqx);
          double maxy = CGAL::abs(pqy);
          double maxz = CGAL::abs(pqz);

          double aprx = CGAL::abs(prx);
          double apry = CGAL::abs(pry);
          double aprz = CGAL::abs(prz);

          double apsx = CGAL::abs(psx);
          double apsy = CGAL::abs(psy);
          double apsz = CGAL::abs(psz);

#ifdef CGAL_USE_SSE2_MAX
          CGAL::Max<double> mmax;
          maxx = mmax(maxx, aqtx, artx, astx);
          maxy = mmax(maxy, aqty, arty, asty);
          maxz = mmax(maxz, aqtz, artz, astz);
#else
          if (maxx < aprx) maxx = aprx;
          if (maxx < apsx) maxx = apsx;

          if (maxy < apry) maxy = apry;
          if (maxy < apsy) maxy = apsy;

          if (maxz < aprz) maxz = aprz;
          if (maxz < apsz) maxz = apsz;
#endif
          double det = CGAL::determinant(pqx, pqy, pqz,
                                         prx, pry, prz,
                                         psx, psy, psz);

          double eps = 4.111024169857068197e-15 * maxx * maxy * maxz;

#ifdef CGAL_USE_SSE2_MAX
          /*
          CGAL::Min<double> mmin;
          double tmp = mmin(maxx, maxy, maxz);
          maxz = mmax(maxx, maxy, maxz);
          maxx = tmp;
          */
          sse2minmax(maxx,maxy,maxz);
          // maxy can contain ANY element
#else
          // Sort maxx < maxy < maxz.
          if (maxx > maxz)
              std::swap(maxx, maxz);
          if (maxy > maxz)
              std::swap(maxy, maxz);
          else if (maxy < maxx)
              std::swap(maxx, maxy);
#endif

          // Protect against underflow in the computation of eps.
          if (maxx < 1e-97) /* cbrt(min_double/eps) */ {
            if (maxx == 0)
              return ZERO;
          }
          // Protect against overflow in the computation of det.
          else if (maxz < 1e102) /* cbrt(max_double [hadamard]/4) */ {
            if (det > eps)  return POSITIVE;
            if (det < -eps) return NEGATIVE;
          }

          CGAL_PROFILER("Periodic_3_orientation_3 semi-static failures");
      }

      return Base::operator()(p,q,r,s,o_p,o_q,o_r,o_s);
  }

  // Computes the epsilon for Periodic_3_orientation_3.
  static double compute_epsilon()
  {
    typedef Static_filter_error F;
    F t1 = F(1, F::ulp()/4);         // First translation
    F det = CGAL::determinant(t1, t1, t1,
                              t1, t1, t1,
                              t1, t1, t1); // Full det
    double err = det.error();
    err += err * 2 * F::ulp(); // Correction due to "eps * maxx * maxy...".
    std::cerr << "*** epsilon for Periodic_3_orientation_3 = " << err
              << std::endl;
    return err;
  }

};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_3_ORIENTATION_3_H
