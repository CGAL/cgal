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

#ifndef CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_4_H
#define CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_4_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Static_filter_error.h>
#include <CGAL/determinant.h>
#include <cmath>

namespace CGAL { namespace internal { namespace Static_filters_predicates {



template < typename K_base >
class Orientation_4
  : public K_base::Orientation_4
{
  typedef typename K_base::Orientation      Orientation;
  typedef typename K_base::Point_4          Point_4;
  typedef typename K_base::Orientation_4    Base;

public:
  using Base::operator();

  Orientation
  operator()(const Point_4 &p, const Point_4 &q,
             const Point_4 &r, const Point_4 &s,
             const Point_4 &t) const
  {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Orientation_4", tmp);

      double p0, p1, p2, p3, q0, q1, q2, q3, r0, r1, r2, r3, s0, s1, s2, s3, t0, t1, t2, t3;

      if (fit_in_double(p.c0(), p0) && fit_in_double(p.c1(), p1) &&
          fit_in_double(p.c2(), p2) && fit_in_double(p.c3(), p3) &&
          fit_in_double(q.c0(), q0) && fit_in_double(q.c1(), q1) &&
          fit_in_double(q.c2(), q2) && fit_in_double(q.c3(), q3) &&
          fit_in_double(r.c0(), r0) && fit_in_double(r.c1(), r1) &&
          fit_in_double(r.c2(), r2) && fit_in_double(r.c3(), r3) &&
          fit_in_double(s.c0(), s0) && fit_in_double(s.c1(), s1) &&
          fit_in_double(s.c2(), s2) && fit_in_double(s.c3(), s3) &&
          fit_in_double(t.c0(), t0) && fit_in_double(t.c1(), t1) &&
          fit_in_double(t.c2(), t2) && fit_in_double(t.c3(), t3)  )
      {
        CGAL_assertion_code(Orientation should_be = Base::operator()(p, q, r, s, t));
double m01;
    m01 = (q0 - p0);
    double m02;
    m02 = (r0 - p0);
    double m03;
    m03 = (s0 - p0);
    double m04;
    m04 = (t0 - p0);
    double m11;
    m11 = (q1 - p1);
    double m12;
    m12 = (r1 - p1);
    double m13;
    m13 = (s1 - p1);
    double m14;
    m14 = (t1 - p1);
    double m21;
    m21 = (q2 - p2);
    double m22;
    m22 = (r2 - p2);
    double m23;
    m23 = (s2 - p2);
    double m24;
    m24 = (t2 - p2);
    double m31;
    m31 = (q3 - p3);
    double m32;
    m32 = (r3 - p3);
    double m33;
    m33 = (s3 - p3);
    double m34;
    m34 = (t3 - p3);
    double det;
    det = ::CGAL::determinant( m01, m02, m03, m04, m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34 );
    int int_tmp_result;
    double eps;
    double max1 = CGAL::abs(m01);
    double am = CGAL::abs(m02);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m03);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m11);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m12);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m13);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m23);
    if( (max1 < am) ) { max1 = am; }


    double max2 = CGAL::abs(m01);
    am = CGAL::abs(m02);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m11);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m12);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m21);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m22);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m23);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m33);
    if( (max2 < am) ) { max2 = am; }


    double max3 = CGAL::abs(m04);
    am = CGAL::abs(m14);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m24);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m34);
    if( (max3 < am) ) { max3 = am; }


    double max4 = CGAL::abs(m11);
    am = CGAL::abs(m12);
    if( (max4 < am) ) { max4 = am; }
    am = CGAL::abs(m21);
    if( (max4 < am) ) { max4 = am; }
    am = CGAL::abs(m22);
    if( (max4 < am) ) { max4 = am; }
    am = CGAL::abs(m31);
    if( (max4 < am) ) { max4 = am; }
    am = CGAL::abs(m32);
    if( (max4 < am) ) { max4 = am; }


    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
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
    if( (lower_bound_1 < 2.89273249588395272840e-74) )
    {
        return Base::operator()(p, q, r, s, t);
    }
    else
    {
        if( (upper_bound_1 > 7.23700557733225900010e+75) )
        {
            return Base::operator()(p, q, r, s, t);
        }
        eps = (3.17768858673611327578e-14 * (((max4 * max2) * max1) * max3));
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
            else
            {
                return Base::operator()(p, q, r, s, t);
            }
        }
    }
}
    return Base::operator()(p, q, r, s, t);

  }

};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_4_H
