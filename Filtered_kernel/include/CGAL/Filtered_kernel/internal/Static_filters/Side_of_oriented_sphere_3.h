// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_INTERNAL_STATIC_FILTERS_SIDE_OF_ORIENTED_SPHERE_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_SIDE_OF_ORIENTED_SPHERE_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Static_filter_error.h>
#include <immintrin.h>

namespace CGAL { namespace internal { namespace Static_filters_predicates {


  inline double hmax4(__m256d v)
{
    __m128d hi = _mm256_extractf128_pd(v, 1);
    __m128d lo = _mm256_castpd256_pd128(v);
    __m128d m  = _mm_max_pd(lo, hi);
    __m128d s  = _mm_shuffle_pd(m, m, 1);
    m = _mm_max_sd(m, s);
    return _mm_cvtsd_f64(m);
}

inline __m256d abs4(__m256d v)
{
    const __m256d mask =
        _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fffffffffffffffLL));
    return _mm256_and_pd(v, mask);
}


template < typename K_base >
class Side_of_oriented_sphere_3
  : public K_base::Side_of_oriented_sphere_3
{
  typedef typename K_base::Oriented_side                Oriented_side;
  typedef typename K_base::Point_3                      Point_3;
  typedef typename K_base::Side_of_oriented_sphere_3    Base;

public:
  Oriented_side
  operator()(const Point_3 &p, const Point_3 &q, const Point_3 &r,
             const Point_3 &s, const Point_3 &t) const
  {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Side_of_oriented_sphere_3", tmp);

      double px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz;

      if (fit_in_double(p.x(), px) && fit_in_double(p.y(), py) &&
          fit_in_double(p.z(), pz) &&
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy) &&
          fit_in_double(q.z(), qz) &&
          fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry) &&
          fit_in_double(r.z(), rz) &&
          fit_in_double(s.x(), sx) && fit_in_double(s.y(), sy) &&
          fit_in_double(s.z(), sz) &&
          fit_in_double(t.x(), tx) && fit_in_double(t.y(), ty) &&
          fit_in_double(t.z(), tz))
      {
          CGAL_BRANCH_PROFILER_BRANCH_1(tmp);

#ifdef CGAL_VECTORIZE
          double X[4] = {px-tx, qx - tx, rx - tx, sx - tx};
          double Y[4] = {py-ty, qy - ty, ry - ty, sy - ty};
          double Z[4] = {pz-tz, qz - tz, rz - tz, sz - tz};

          double pt2 = CGAL_NTS square(X[0]) + CGAL_NTS square(Y[0])
                     + CGAL_NTS square(Z[0]);
          double qt2 = CGAL_NTS square(X[1]) + CGAL_NTS square(Y[1])
                     + CGAL_NTS square(Z[1]);
          double rt2 = CGAL_NTS square(X[2]) + CGAL_NTS square(Y[2])
                     + CGAL_NTS square(Z[2]);
          double st2 = CGAL_NTS square(X[3]) + CGAL_NTS square(Y[3])
                     + CGAL_NTS square(Z[3]);

#else
          double ptx = px - tx;
          double pty = py - ty;
          double ptz = pz - tz;

          double qtx = qx - tx;
          double qty = qy - ty;
          double qtz = qz - tz;

          double rtx = rx - tx;
          double rty = ry - ty;
          double rtz = rz - tz;

          double stx = sx - tx;
          double sty = sy - ty;
          double stz = sz - tz;


          double pt2 = CGAL_NTS square(ptx) + CGAL_NTS square(pty)
                     + CGAL_NTS square(ptz);

          double qt2 = CGAL_NTS square(qtx) + CGAL_NTS square(qty)
                     + CGAL_NTS square(qtz);
          double rt2 = CGAL_NTS square(rtx) + CGAL_NTS square(rty)
                     + CGAL_NTS square(rtz);
          double st2 = CGAL_NTS square(stx) + CGAL_NTS square(sty)
                     + CGAL_NTS square(stz);
#endif
          // Compute the semi-static bound.

#if CGAL_STD_ABS
double maxx = (std::max)({std::abs(ptx), std::abs(qtx), std::abs(rtx), std::abs(stx)});
double maxy = (std::max)({std::abs(pty), std::abs(qty), std::abs(rty), std::abs(sty)});
double maxz = (std::max)({std::abs(ptz), std::abs(qtz), std::abs(rtz), std::abs(stz)});
#elif CGAL_STD_FABS
double maxx = (std::max)({std::fabs(ptx), std::fabs(qtx), std::fabs(rtx), std::fabs(stx)});
double maxy = (std::max)({std::fabs(pty), std::fabs(qty), std::fabs(rty), std::fabs(sty)});
double maxz = (std::max)({std::fabs(ptz), std::fabs(qtz), std::fabs(rtz), std::fabs(stz)});

#elif CGAL_VECTORIZE
double maxx = hmax4(abs4(_mm256_loadu_pd(X)));
double maxy = hmax4(abs4(_mm256_loadu_pd(Y)));
double maxz = hmax4(abs4(_mm256_loadu_pd(Z)));
#elif CGAL_CGAL_ABS
double maxx = (std::max)({CGAL::abs(ptx), CGAL::abs(qtx), CGAL::abs(rtx), CGAL::abs(stx)});
double maxy = (std::max)({CGAL::abs(pty), CGAL::abs(qty), CGAL::abs(rty), CGAL::abs(sty)});
double maxz = (std::max)({CGAL::abs(ptz), CGAL::abs(qtz), CGAL::abs(rtz), CGAL::abs(stz)});
#else
          double maxx = CGAL::abs(ptx);
          double maxy = CGAL::abs(pty);
          double maxz = CGAL::abs(ptz);

          double aqtx = CGAL::abs(qtx);
          double artx = CGAL::abs(rtx);
          double astx = CGAL::abs(stx);

          double aqty = CGAL::abs(qty);
          double arty = CGAL::abs(rty);
          double asty = CGAL::abs(sty);

          double aqtz = CGAL::abs(qtz);
          double artz = CGAL::abs(rtz);
          double astz = CGAL::abs(stz);


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
#ifdef CGAL_VECTORIZE
          double det = CGAL::determinant(X[0],Y[0],Z[0],pt2,
                                         X[2],Y[2],Z[2],rt2,
                                         X[1],Y[1],Z[1],qt2,
                                         X[3],Y[3],Z[3],st2);
#else
          double det = CGAL::determinant(ptx,pty,ptz,pt2,
                                         rtx,rty,rtz,rt2,
                                         qtx,qty,qtz,qt2,
                                         stx,sty,stz,st2);
#endif
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

          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
      }
      return Base::operator()(p, q, r, s, t);
  }

  // Computes the epsilon for Side_of_oriented_sphere_3.
  static double compute_epsilon()
  {
    typedef internal::Static_filter_error F;
    F t1 = F(1,F::ulp()/2);         // First translation
    F sq = t1*t1+t1*t1+t1*t1; // squares
    F det = CGAL::determinant(t1, t1, t1, sq,
                              t1, t1, t1, sq,
                              t1, t1, t1, sq,
                              t1, t1, t1, sq); // Full det
    double err = det.error();
    err += err * 3 * F::ulp(); // Correction due to "eps * maxx * ...".

    std::cerr << "*** epsilon for Side_of_oriented_sphere_3 = " << err << std::endl;
    return err;
  }
};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_SIDE_OF_ORIENTED_SPHERE_3_H
