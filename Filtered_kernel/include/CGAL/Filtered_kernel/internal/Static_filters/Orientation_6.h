// Copyright (c) 2025  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_6_H
#define CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_6_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Static_filter_error.h>
#include <CGAL/determinant.h>
#include <cmath>

namespace CGAL { namespace internal { namespace Static_filters_predicates {



template < typename K_base >
class Orientation_6
  : public K_base::Orientation_6
{
  typedef typename K_base::Orientation      Orientation;
  typedef typename K_base::Point_6          Point_6;
  typedef typename K_base::Orientation_6    Base;

public:
  using Base::operator();

  Orientation
  operator()(const Point_6 &p, const Point_6 &q,
             const Point_6 &r, const Point_6 &s,
             const Point_6 &t, const Point_6 &u, const Point_6 &v) const
  {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Orientation_6", tmp);

      double p0, p1, p2, p3, p4, p5, q0, q1, q2, q3, q4, q5, r0, r1, r2, r3, r4, r5, s0, s1, s2, s3, s4, s5, t0, t1, t2, t3, t4, t5, u0, u1, u2, u3, u4, u5, v0, v1, v2, v3, v4, v5;
      if (fit_in_double(p.c0(), p0) && fit_in_double(p.c1(), p1) &&
          fit_in_double(p.c2(), p2) && fit_in_double(p.c3(), p3) &&
          fit_in_double(p.c4(), p4) && fit_in_double(p.c5(), p5) &&

          fit_in_double(q.c0(), q0) && fit_in_double(q.c1(), q1) &&
          fit_in_double(q.c2(), q2) && fit_in_double(q.c3(), q3) &&
          fit_in_double(q.c4(), q4) && fit_in_double(q.c5(), q5) &&

          fit_in_double(r.c0(), r0) && fit_in_double(r.c1(), r1) &&
          fit_in_double(r.c2(), r2) && fit_in_double(r.c3(), r3) &&
          fit_in_double(r.c4(), r4) && fit_in_double(r.c5(), r5) &&

          fit_in_double(s.c0(), s0) && fit_in_double(s.c1(), s1) &&
          fit_in_double(s.c2(), s2) && fit_in_double(s.c3(), s3) &&
          fit_in_double(s.c4(), s4) && fit_in_double(s.c5(), s5) &&

          fit_in_double(t.c0(), t0) && fit_in_double(t.c1(), t1) &&
          fit_in_double(t.c2(), t2) && fit_in_double(t.c3(), t3) &&
          fit_in_double(t.c4(), t4) && fit_in_double(t.c5(), t5) &&

          fit_in_double(u.c0(), u0) && fit_in_double(u.c1(), u1) &&
          fit_in_double(u.c2(), u2) && fit_in_double(u.c3(), u3) &&
          fit_in_double(u.c4(), u4) && fit_in_double(u.c5(), u5) &&

          fit_in_double(v.c0(), v0) && fit_in_double(v.c1(), v1) &&
          fit_in_double(v.c2(), v2) && fit_in_double(v.c3(), v3) &&
          fit_in_double(v.c4(), v4) && fit_in_double(v.c5(), v5))
      {
        CGAL_assertion_code(Orientation should_be = Base::operator()(p, q, r, s, t, u, v));
        double m01;
    m01 = (q0 - p0);
    double m02;
    m02 = (r0 - p0);
    double m03;
    m03 = (s0 - p0);
    double m04;
    m04 = (t0 - p0);
    double m05;
    m05 = (u0 - p0);
    double m06;
    m06 = (v0 - p0);
    double m11;
    m11 = (q1 - p1);
    double m12;
    m12 = (r1 - p1);
    double m13;
    m13 = (s1 - p1);
    double m14;
    m14 = (t1 - p1);
    double m15;
    m15 = (u1 - p1);
    double m16;
    m16 = (v1 - p1);
    double m21;
    m21 = (q2 - p2);
    double m22;
    m22 = (r2 - p2);
    double m23;
    m23 = (s2 - p2);
    double m24;
    m24 = (t2 - p2);
    double m25;
    m25 = (u2 - p2);
    double m26;
    m26 = (v2 - p2);
    double m31;
    m31 = (q3 - p3);
    double m32;
    m32 = (r3 - p3);
    double m33;
    m33 = (s3 - p3);
    double m34;
    m34 = (t3 - p3);
    double m35;
    m35 = (u3 - p3);
    double m36;
    m36 = (v3 - p3);
    double m41;
    m41 = (q4 - p4);
    double m42;
    m42 = (r4 - p4);
    double m43;
    m43 = (s4 - p4);
    double m44;
    m44 = (t4 - p4);
    double m45;
    m45 = (u4 - p4);
    double m46;
    m46 = (v4 - p4);
    double m51;
    m51 = (q5 - p5);
    double m52;
    m52 = (r5 - p5);
    double m53;
    m53 = (s5 - p5);
    double m54;
    m54 = (t5 - p5);
    double m55;
    m55 = (u5 - p5);
    double m56;
    m56 = (v5 - p5);
    double det;
    det = ::CGAL::determinant( m01, m02, m03, m04, m05, m06, m11, m12, m13, m14, m15, m16, m21, m22, m23, m24, m25, m26, m31, m32, m33, m34, m35, m36, m41, m42, m43, m44, m45, m46, m51, m52, m53, m54, m55, m56 );
    double eps;
    double max1 = CGAL::abs(m01);
    double am = CGAL::abs(m02);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m11);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m12);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m21);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m22);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m31);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m32);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m41);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m42);
    if( (max1 < am) ) { max1 = am; }


    double max2 = CGAL::abs(m03);
    am = CGAL::abs(m11);

    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m12);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m13);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m21);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m22);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m23);

    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m31);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m32);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m33);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m41);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m42);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m43);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m51);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m52);
    if( (max2 < am) ) { max2 = am; }


    double max3 = CGAL::abs(m04);
    am = CGAL::abs(m14);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m24);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m34);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m44);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m54);
    if( (max3 < am) ) { max3 = am; }

    double max4 = CGAL::abs(m05);
    am = CGAL::abs(m15);
    if( (max4 < am) ) { max4 = am; }
    am = CGAL::abs(m25);
    if( (max4 < am) ) { max3 = am; }
    am = CGAL::abs(m35);
    if( (max4 < am) ) { max4= am; }
    am = CGAL::abs(m45);
    if( (max4 < am) ) { max4 = am; }
    am = CGAL::abs(m55);
    if( (max4 < am) ) { max4 = am; }

    double max5 = CGAL::abs(m06);
    am = CGAL::abs(m16);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m26);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m36);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m46);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m56);
    if( (max5 < am) ) { max5 = am; }


    double max6 = CGAL::abs(m13);
    am = CGAL::abs(m21);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m22);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m23);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m31);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m32);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m33);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m41);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m42);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m43);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m51);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m52);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m53);
    am = CGAL::abs(m53);
    if( (max6 < am) ) { max6 = am; }


    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max6;
    upper_bound_1 = max6;
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (lower_bound_1 < 4.82472686053427214432e-50) )
    {
        return Base::operator()(p, q, r, s, t, u, v);
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355511223e+50) )
        {
            return Base::operator()(p, q, r, s, t, u, v);
        }
        eps = (1.76403842114300859158e-12 * (((((max1 * max2) * max6) * max3) * max4) * max5));
        if( (det > eps) )
        {
            CGAL_assertion(should_be == POSITIVE);
            return POSITIVE;
        }
        else
        {
            if( (det < -eps) )
            {
                CGAL_assertion(should_be == NEGATIVE);
                return NEGATIVE;
            }
        }
    }
}
    return Base::operator()(p, q, r, s, t, u, v);
  }
};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_6_H
