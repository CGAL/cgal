// Copyright (c) 2001,2004,2008-2009   INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_3_SIDE_OF_ORIENTED_SPHERE_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_3_SIDE_OF_ORIENTED_SPHERE_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <CGAL/internal/Static_filters/tools.h>

#include <CGAL/Periodic_3_offset_3.h>

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {

template < class K, class Side_of_oriented_sphere_3_base >
class Periodic_3_side_of_oriented_sphere_3
  : public Side_of_oriented_sphere_3_base
{
  typedef Side_of_oriented_sphere_3_base                Base;

public:
  typedef K                                             Kernel;

  typedef typename K::FT                                FT;
  typedef typename K::Point_3                           Point_3;
  typedef typename K::Iso_cuboid_3                      Iso_cuboid_3;
  typedef typename K::Periodic_3_offset_3               Offset;

public:
  typedef typename Base::result_type  result_type;

private:
  const Iso_cuboid_3 * _dom;

public:
  Periodic_3_side_of_oriented_sphere_3(const Iso_cuboid_3 * dom,
                                       const Side_of_oriented_sphere_3_base& sop3)
      : Base(sop3), _dom(dom)
  { }

  Oriented_side
  operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r,
             const Point_3 &s, const Point_3 &t) const
  {
      CGAL_PROFILER("Periodic_3_side_of_oriented_sphere_3 calls");

      Get_approx<Point_3> get_approx; // Identity functor for all points
                                      // but lazy points.
      double px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz;
      init_double(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, (FT*)(0));

      if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
          fit_in_double(get_approx(p).z(), pz) &&
          fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
          fit_in_double(get_approx(q).z(), qz) &&
          fit_in_double(get_approx(r).x(), rx) && fit_in_double(get_approx(r).y(), ry) &&
          fit_in_double(get_approx(r).z(), rz) &&
          fit_in_double(get_approx(s).x(), sx) && fit_in_double(get_approx(s).y(), sy) &&
          fit_in_double(get_approx(s).z(), sz) &&
          fit_in_double(get_approx(t).x(), tx) && fit_in_double(get_approx(t).y(), ty) &&
          fit_in_double(get_approx(t).z(), tz))
      {
          CGAL_PROFILER("Periodic_3_side_of_oriented_sphere_3 semi-static attempts");

          double ptx = px - tx;
          double pty = py - ty;
          double ptz = pz - tz;
          double pt2 = CGAL_NTS square(ptx) + CGAL_NTS square(pty)
                       + CGAL_NTS square(ptz);
          double qtx = qx - tx;
          double qty = qy - ty;
          double qtz = qz - tz;
          double qt2 = CGAL_NTS square(qtx) + CGAL_NTS square(qty)
                       + CGAL_NTS square(qtz);
          double rtx = rx - tx;
          double rty = ry - ty;
          double rtz = rz - tz;
          double rt2 = CGAL_NTS square(rtx) + CGAL_NTS square(rty)
                       + CGAL_NTS square(rtz);
          double stx = sx - tx;
          double sty = sy - ty;
          double stz = sz - tz;
          double st2 = CGAL_NTS square(stx) + CGAL_NTS square(sty)
                       + CGAL_NTS square(stz);

          // Compute the semi-static bound.
          double maxx = CGAL::abs(ptx);
          double maxy = CGAL::abs(pty);
          double maxz = CGAL::abs(ptz);

          double aqtx = CGAL::abs(qtx);
          double aqty = CGAL::abs(qty);
          double aqtz = CGAL::abs(qtz);

          double artx = CGAL::abs(rtx);
          double arty = CGAL::abs(rty);
          double artz = CGAL::abs(rtz);

          double astx = CGAL::abs(stx);
          double asty = CGAL::abs(sty);
          double astz = CGAL::abs(stz);

#ifdef CGAL_USE_SSE2_MAX
          CGAL::Max<double> mmax;
          maxx = mmax(maxx, aqtx, artx, astx);
          maxy = mmax(maxy, aqty, arty, asty);
          maxz = mmax(maxz, aqtz, artz, astz);
#else
          if (maxx < aqtx) maxx = aqtx;
          if (maxx < artx) maxx = artx;
          if (maxx < astx) maxx = astx;

          if (maxy < aqty) maxy = aqty;
          if (maxy < arty) maxy = arty;
          if (maxy < asty) maxy = asty;

          if (maxz < aqtz) maxz = aqtz;
          if (maxz < artz) maxz = artz;
          if (maxz < astz) maxz = astz;
#endif

          double eps = 1.2466136531027298e-13 * maxx * maxy * maxz;

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
          double det = CGAL::determinant(ptx,pty,ptz,pt2,
                                         rtx,rty,rtz,rt2,
                                         qtx,qty,qtz,qt2,
                                         stx,sty,stz,st2);

          // Protect against underflow in the computation of eps.
          if (maxx < 1e-58) /* sqrt^5(min_double/eps) */ {
            if (maxx == 0)
              return ON_ORIENTED_BOUNDARY;
          }
          // Protect against overflow in the computation of det.
          else if (maxz < 1e61) /* sqrt^5(max_double/4 [hadamard]) */ {
            eps *= (maxz * maxz);
            if (det > eps)  return ON_POSITIVE_SIDE;
            if (det < -eps) return ON_NEGATIVE_SIDE;
          }

          CGAL_PROFILER("Periodic_3_side_of_oriented_sphere_3 semi-static failures");
      }
      return Base::operator()(p, q, r, s, t);
  }

  Oriented_side
  operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r,
             const Point_3 &s, const Point_3 &t,
             const Offset &o_p, const Offset &o_q, const Offset &o_r,
             const Offset &o_s, const Offset &o_t) const
  {
    CGAL_PROFILER("Periodic_3_side_of_oriented_sphere_3 calls");

    Get_approx<Point_3> get_approx; // Identity functor for all points
                                    // but lazy points.

    double px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz;
    double domxmax, domxmin, domymax, domymin, domzmax, domzmin;
    init_double(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, (FT*)(0));
    init_double(domxmax, domxmin, domymax, domymin, domzmax, domzmin, (FT*)(0));
    int otx = o_t.x();
    int oty = o_t.y();
    int otz = o_t.z();

    if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
        fit_in_double(get_approx(p).z(), pz) &&
        fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
        fit_in_double(get_approx(q).z(), qz) &&
        fit_in_double(get_approx(r).x(), rx) && fit_in_double(get_approx(r).y(), ry) &&
        fit_in_double(get_approx(r).z(), rz) &&
        fit_in_double(get_approx(s).x(), sx) && fit_in_double(get_approx(s).y(), sy) &&
        fit_in_double(get_approx(s).z(), sz) &&
        fit_in_double(get_approx(t).x(), tx) && fit_in_double(get_approx(t).y(), ty) &&
        fit_in_double(get_approx(t).z(), tz) &&
        fit_in_double(_dom->xmax(), domxmax) &&
        fit_in_double(_dom->xmin(), domxmin) &&
        fit_in_double(_dom->ymax(), domymax) &&
        fit_in_double(_dom->ymin(), domymin) &&
        fit_in_double(_dom->zmax(), domzmax) &&
        fit_in_double(_dom->zmin(), domzmin))
    {
      CGAL_PROFILER("Periodic_3_side_of_oriented_sphere_3 semi-static attempts");

      double domx = domxmax - domxmin;
      double domy = domymax - domymin;
      double domz = domzmax - domzmin;

      double ptx = px - tx + domx * (o_p.x() - otx);
      double pty = py - ty + domy * (o_p.y() - oty);
      double ptz = pz - tz + domz * (o_p.z() - otz);
      double pt2 = CGAL_NTS square(ptx) + CGAL_NTS square(pty)
            + CGAL_NTS square(ptz);
      double qtx = qx - tx + domx * (o_q.x() - otx);
      double qty = qy - ty + domy * (o_q.y() - oty);
      double qtz = qz - tz + domz * (o_q.z() - otz);
      double qt2 = CGAL_NTS square(qtx) + CGAL_NTS square(qty)
            + CGAL_NTS square(qtz);
      double rtx = rx - tx + domx * (o_r.x() - otx);
      double rty = ry - ty + domy * (o_r.y() - oty);
      double rtz = rz - tz + domz * (o_r.z() - otz);
      double rt2 = CGAL_NTS square(rtx) + CGAL_NTS square(rty)
            + CGAL_NTS square(rtz);
      double stx = sx - tx + domx * (o_s.x() - otx);
      double sty = sy - ty + domy * (o_s.y() - oty);
      double stz = sz - tz + domz * (o_s.z() - otz);
      double st2 = CGAL_NTS square(stx) + CGAL_NTS square(sty)
            + CGAL_NTS square(stz);

      // Compute the semi-static bound.
      double maxx = CGAL::abs(ptx);
      double maxy = CGAL::abs(pty);
      double maxz = CGAL::abs(ptz);

      double aqtx = CGAL::abs(qtx);
      double aqty = CGAL::abs(qty);
      double aqtz = CGAL::abs(qtz);

      double artx = CGAL::abs(rtx);
      double arty = CGAL::abs(rty);
      double artz = CGAL::abs(rtz);

      double astx = CGAL::abs(stx);
      double asty = CGAL::abs(sty);
      double astz = CGAL::abs(stz);

#ifdef CGAL_USE_SSE2_MAX
          CGAL::Max<double> mmax;
          maxx = mmax(maxx, aqtx, artx, astx);
          maxy = mmax(maxy, aqty, arty, asty);
          maxz = mmax(maxz, aqtz, artz, astz);
#else
      if (maxx < aqtx) maxx = aqtx;
      if (maxx < artx) maxx = artx;
      if (maxx < astx) maxx = astx;

      if (maxy < aqty) maxy = aqty;
      if (maxy < arty) maxy = arty;
      if (maxy < asty) maxy = asty;

      if (maxz < aqtz) maxz = aqtz;
      if (maxz < artz) maxz = artz;
      if (maxz < astz) maxz = astz;
#endif

      double eps = 1.0466759304746772485e-13 * maxx * maxy * maxz;

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
      double det = CGAL::determinant(ptx,pty,ptz,pt2,
                                     rtx,rty,rtz,rt2,
                                     qtx,qty,qtz,qt2,
                                     stx,sty,stz,st2);

      // Protect against underflow in the computation of eps.
      if (maxx < 1e-58) /* sqrt^5(min_double/eps) */ {
        if (maxx == 0)
          return ON_ORIENTED_BOUNDARY;
      }
      // Protect against overflow in the computation of det.
      else if (maxz < 1e61) /* sqrt^5(max_double/4 [hadamard]) */ {
        eps *= (maxz * maxz);
        if (det > eps)  return ON_POSITIVE_SIDE;
        if (det < -eps) return ON_NEGATIVE_SIDE;
      }

      CGAL_PROFILER("Periodic_3_side_of_oriented_sphere_3 semi-static failures");
    }
    return Base::operator()(p, q, r, s, t, o_p, o_q, o_r, o_s, o_t);
  }

  // Computes the epsilon for Periodic_3_side_of_oriented_sphere_3.
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

    std::cerr << "*** epsilon for Periodic_3_side_of_oriented_sphere_3 = "
              << err << std::endl;
    return err;
  }
};

} // namespace Static_filters_predicates

} // namespace internal

} // namespace CGAL

#endif // CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_3_SIDE_OF_ORIENTED_SPHERE_3_H
