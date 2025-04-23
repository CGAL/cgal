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

#ifndef CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_5_H
#define CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_5_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Static_filter_error.h>
#include <cmath>

namespace CGAL { namespace internal { namespace Static_filters_predicates {



template < typename K_base >
class Orientation_5
  : public K_base::Orientation_5
{
  typedef typename K_base::Orientation      Orientation;
  typedef typename K_base::Point_5          Point_5;
  typedef typename K_base::Orientation_5    Base;

public:
  using Base::operator();

  Orientation
  operator()(const Point_5 &p, const Point_5 &q,
             const Point_5 &r, const Point_5 &s,
             const Point_5 &t, const Point_5 &u) const
  {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Orientation_5", tmp);

      double p0, p1, p2, p3, p4, q0, q1, q2, q3, q4, r0, r1, r2, r3, r4, s0, s1, s2, s3, s4, t0, t1, t2, t3, t4, u0, u1, u2, u3, u4;

      if (fit_in_double(p.c0(), p0) && fit_in_double(p.c1(), p1) &&
          fit_in_double(p.c2(), p2) && fit_in_double(p.c3(), p3) &&
          fit_in_double(p.c4(), p4) &&
          fit_in_double(q.c0(), q0) && fit_in_double(q.c1(), q1) &&
          fit_in_double(q.c2(), q2) && fit_in_double(q.c3(), q3) &&
          fit_in_double(q.c4(), q4) &&
          fit_in_double(r.c0(), r0) && fit_in_double(r.c1(), r1) &&
          fit_in_double(r.c2(), r2) && fit_in_double(r.c3(), r3) &&
          fit_in_double(r.c4(), r4) &&
          fit_in_double(s.c0(), s0) && fit_in_double(s.c1(), s1) &&
          fit_in_double(s.c2(), s2) && fit_in_double(s.c3(), s3) &&
          fit_in_double(s.c4(), s4) &&
          fit_in_double(t.c0(), t0) && fit_in_double(t.c1(), t1) &&
          fit_in_double(t.c2(), t2) && fit_in_double(t.c3(), t3) &&
          fit_in_double(t.c4(), t4) &&
          fit_in_double(u.c0(), u0) && fit_in_double(u.c1(), u1) &&
          fit_in_double(u.c2(), u2) && fit_in_double(u.c3(), u3) &&
          fit_in_double(u.c4(), u4))
      {
        CGAL_assertion_code(Orientation should_be = Base::operator()(p, q, r, s, t, u));
        double det;
    double determinant_with_1_in_row_0_return_value;
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
    double m12;
    m12 = (r0 - q0);
    double m13;
    m13 = (s0 - q0);
    double m14;
    m14 = (t0 - q0);
    double m15;
    m15 = (u0 - q0);
    double m23;
    m23 = (s0 - r0);
    double m24;
    m24 = (t0 - r0);
    double m25;
    m25 = (u0 - r0);
    double m34;
    m34 = (t0 - s0);
    double m35;
    m35 = (u0 - s0);
    double m45;
    m45 = (u0 - t0);
    double m012;
    m012 = (((m01 * r1) - (m02 * q1)) + (m12 * p1));
    double m013;
    m013 = (((m01 * s1) - (m03 * q1)) + (m13 * p1));
    double m014;
    m014 = (((m01 * t1) - (m04 * q1)) + (m14 * p1));
    double m015;
    m015 = (((m01 * u1) - (m05 * q1)) + (m15 * p1));
    double m023;
    m023 = (((m02 * s1) - (m03 * r1)) + (m23 * p1));
    double m024;
    m024 = (((m02 * t1) - (m04 * r1)) + (m24 * p1));
    double m025;
    m025 = (((m02 * u1) - (m05 * r1)) + (m25 * p1));
    double m034;
    m034 = (((m03 * t1) - (m04 * s1)) + (m34 * p1));
    double m035;
    m035 = (((m03 * u1) - (m05 * s1)) + (m35 * p1));
    double m045;
    m045 = (((m04 * u1) - (m05 * t1)) + (m45 * p1));
    double m123;
    m123 = (((m12 * s1) - (m13 * r1)) + (m23 * q1));
    double m124;
    m124 = (((m12 * t1) - (m14 * r1)) + (m24 * q1));
    double m125;
    m125 = (((m12 * u1) - (m15 * r1)) + (m25 * q1));
    double m134;
    m134 = (((m13 * t1) - (m14 * s1)) + (m34 * q1));
    double m135;
    m135 = (((m13 * u1) - (m15 * s1)) + (m35 * q1));
    double m145;
    m145 = (((m14 * u1) - (m15 * t1)) + (m45 * q1));
    double m234;
    m234 = (((m23 * t1) - (m24 * s1)) + (m34 * r1));
    double m235;
    m235 = (((m23 * u1) - (m25 * s1)) + (m35 * r1));
    double m245;
    m245 = (((m24 * u1) - (m25 * t1)) + (m45 * r1));
    double m345;
    m345 = (((m34 * u1) - (m35 * t1)) + (m45 * s1));
    double m0123;
    m0123 = ((((m012 * s2) - (m013 * r2)) + (m023 * q2)) - (m123 * p2));
    double m0124;
    m0124 = ((((m012 * t2) - (m014 * r2)) + (m024 * q2)) - (m124 * p2));
    double m0125;
    m0125 = ((((m012 * u2) - (m015 * r2)) + (m025 * q2)) - (m125 * p2));
    double m0134;
    m0134 = ((((m013 * t2) - (m014 * s2)) + (m034 * q2)) - (m134 * p2));
    double m0135;
    m0135 = ((((m013 * u2) - (m015 * s2)) + (m035 * q2)) - (m135 * p2));
    double m0145;
    m0145 = ((((m014 * u2) - (m015 * t2)) + (m045 * q2)) - (m145 * p2));
    double m0234;
    m0234 = ((((m023 * t2) - (m024 * s2)) + (m034 * r2)) - (m234 * p2));
    double m0235;
    m0235 = ((((m023 * u2) - (m025 * s2)) + (m035 * r2)) - (m235 * p2));
    double m0245;
    m0245 = ((((m024 * u2) - (m025 * t2)) + (m045 * r2)) - (m245 * p2));
    double m0345;
    m0345 = ((((m034 * u2) - (m035 * t2)) + (m045 * s2)) - (m345 * p2));
    double m1234;
    m1234 = ((((m123 * t2) - (m124 * s2)) + (m134 * r2)) - (m234 * q2));
    double m1235;
    m1235 = ((((m123 * u2) - (m125 * s2)) + (m135 * r2)) - (m235 * q2));
    double m1245;
    m1245 = ((((m124 * u2) - (m125 * t2)) + (m145 * r2)) - (m245 * q2));
    double m1345;
    m1345 = ((((m134 * u2) - (m135 * t2)) + (m145 * s2)) - (m345 * q2));
    double m2345;
    m2345 = ((((m234 * u2) - (m235 * t2)) + (m245 * s2)) - (m345 * r2));
    double m01234;
    m01234 = (((((m0123 * t3) - (m0124 * s3)) + (m0134 * r3)) - (m0234 * q3)) + (m1234 * p3));
    double m01235;
    m01235 = (((((m0123 * u3) - (m0125 * s3)) + (m0135 * r3)) - (m0235 * q3)) + (m1235 * p3));
    double m01245;
    m01245 = (((((m0124 * u3) - (m0125 * t3)) + (m0145 * r3)) - (m0245 * q3)) + (m1245 * p3));
    double m01345;
    m01345 = (((((m0134 * u3) - (m0135 * t3)) + (m0145 * s3)) - (m0345 * q3)) + (m1345 * p3));
    double m02345;
    m02345 = (((((m0234 * u3) - (m0235 * t3)) + (m0245 * s3)) - (m0345 * r3)) + (m2345 * p3));
    double m12345;
    m12345 = (((((m1234 * u3) - (m1235 * t3)) + (m1245 * s3)) - (m1345 * r3)) + (m2345 * q3));
    double m012345;
    m012345 = ((((((m01234 * u4) - (m01235 * t4)) + (m01245 * s4)) - (m01345 * r4)) + (m02345 * q4)) - (m12345 * p4));
    determinant_with_1_in_row_0_return_value = m012345;
    det = determinant_with_1_in_row_0_return_value;
    double eps;
    double max1 = CGAL::abs(p1);
    if( (max1 < CGAL::abs(q1)) )
    {
        max1 = CGAL::abs(q1);
    }
    if( (max1 < CGAL::abs(r1)) )
    {
        max1 = CGAL::abs(r1);
    }
    if( (max1 < CGAL::abs(s1)) )
    {
        max1 = CGAL::abs(s1);
    }
    if( (max1 < CGAL::abs(t1)) )
    {
        max1 = CGAL::abs(t1);
    }
    if( (max1 < CGAL::abs(u1)) )
    {
        max1 = CGAL::abs(u1);
    }
    double max2 = CGAL::abs(p2);
    if( (max2 < CGAL::abs(q2)) )
    {
        max2 = CGAL::abs(q2);
    }
    if( (max2 < CGAL::abs(r2)) )
    {
        max2 = CGAL::abs(r2);
    }
    if( (max2 < CGAL::abs(s2)) )
    {
        max2 = CGAL::abs(s2);
    }
    if( (max2 < CGAL::abs(t2)) )
    {
        max2 = CGAL::abs(t2);
    }
    if( (max2 < CGAL::abs(u2)) )
    {
        max2 = CGAL::abs(u2);
    }
    double max3 = CGAL::abs(p3);
    if( (max3 < CGAL::abs(q3)) )
    {
        max3 = CGAL::abs(q3);
    }
    if( (max3 < CGAL::abs(r3)) )
    {
        max3 = CGAL::abs(r3);
    }
    if( (max3 < CGAL::abs(s3)) )
    {
        max3 = CGAL::abs(s3);
    }
    if( (max3 < CGAL::abs(t3)) )
    {
        max3 = CGAL::abs(t3);
    }
    if( (max3 < CGAL::abs(u3)) )
    {
        max3 = CGAL::abs(u3);
    }
    double max4 = CGAL::abs(p4);
    if( (max4 < CGAL::abs(q4)) )
    {
        max4 = CGAL::abs(q4);
    }
    if( (max4 < CGAL::abs(r4)) )
    {
        max4 = CGAL::abs(r4);
    }
    if( (max4 < CGAL::abs(s4)) )
    {
        max4 = CGAL::abs(s4);
    }
    if( (max4 < CGAL::abs(t4)) )
    {
        max4 = CGAL::abs(t4);
    }
    if( (max4 < CGAL::abs(u4)) )
    {
        max4 = CGAL::abs(u4);
    }
    double max5 = CGAL::abs(m01);
    if( (max5 < CGAL::abs(m05)) )
    {
        max5 = CGAL::abs(m05);
    }
    if( (max5 < CGAL::abs(m04)) )
    {
        max5 = CGAL::abs(m04);
    }
    if( (max5 < CGAL::abs(m03)) )
    {
        max5 = CGAL::abs(m03);
    }
    if( (max5 < CGAL::abs(m02)) )
    {
        max5 = CGAL::abs(m02);
    }
    if( (max5 < CGAL::abs(m15)) )
    {
        max5 = CGAL::abs(m15);
    }
    if( (max5 < CGAL::abs(m14)) )
    {
        max5 = CGAL::abs(m14);
    }
    if( (max5 < CGAL::abs(m13)) )
    {
        max5 = CGAL::abs(m13);
    }
    if( (max5 < CGAL::abs(m12)) )
    {
        max5 = CGAL::abs(m12);
    }
    if( (max5 < CGAL::abs(m34)) )
    {
        max5 = CGAL::abs(m34);
    }
    if( (max5 < CGAL::abs(m25)) )
    {
        max5 = CGAL::abs(m25);
    }
    if( (max5 < CGAL::abs(m24)) )
    {
        max5 = CGAL::abs(m24);
    }
    if( (max5 < CGAL::abs(m23)) )
    {
        max5 = CGAL::abs(m23);
    }
    if( (max5 < CGAL::abs(m45)) )
    {
        max5 = CGAL::abs(m45);
    }
    if( (max5 < CGAL::abs(m35)) )
    {
        max5 = CGAL::abs(m35);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max5;
    upper_bound_1 = max5;
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
    if( (lower_bound_1 < 8.19482853969781542511e-60) )
    {
        return Base::operator()(p, q, r, s, t, u);
    }
    else
    {
        if( (upper_bound_1 > 3.21387608851797912384e+60) )
        {
            return Base::operator()(p, q, r, s, t, u);
        }
        eps = (6.02067348555779570000e-13 * ((((max5 * max1) * max2) * max3) * max4));
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
      return Base::operator()(p, q, r, s, t, u);
  }


};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_5_H
