// Copyright (c) 20025  GeometryFactory (France).
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

    // return Base::operator()(p, q, r, s, t, u);
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
        Orientation should_be = Base::operator()(p, q, r, s, t, u);
        double pq0;
        pq0 = (q0 - p0);
        double pq1;
        pq1 = (q1 - p1);
        double pq2;
        pq2 = (q2 - p2);
        double pq3;
        pq3 = (q3 - p3);
        double pq4;
        pq4 = (q4 - p4);
        double pr0;
        pr0 = (r0 - p0);
        double pr1;
        pr1 = (r1 - p1);
        double pr2;
        pr2 = (r2 - p2);
        double pr3;
        pr3 = (r3 - p3);
        double pr4;
        pr4 = (r4 - p4);
        double ps0;
        ps0 = (s0 - p0);
        double ps1;
        ps1 = (s1 - p1);
        double ps2;
        ps2 = (s2 - p2);
        double ps3;
        ps3 = (s3 - p3);
        double ps4;
        ps4 = (s4 - p4);
        double pt0;
        pt0 = (t0 - p0);
        double pt1;
        pt1 = (t1 - p1);
        double pt2;
        pt2 = (t2 - p2);
        double pt3;
        pt3 = (t3 - p3);
        double pt4;
        pt4 = (t4 - p4);
        double pu0;
        pu0 = (t0 - p0);
        double pu1;
        pu1 = (u1 - p1);
        double pu2;
        pu2 = (u2 - p2);
        double pu3;
        pu3 = (u3 - p3);
        double pu4;
        pu4 = (u4 - p4);
        double det;
        det = determinant( pq0, pq1, pq2, pq3, pq4, pr0, pr1, pr2, pr3, pr4, ps0, ps1, ps2, ps3, ps4, pt0, pt1, pt2, pt3, pt4, pu0, pu1, pu2, pu3, pu4 );
        Orientation result = ZERO;
        double eps;
        double max1 = CGAL::abs(pq0);
        if( (max1 < CGAL::abs(pr0)) )
        {
            max1 = CGAL::abs(pr0);
        }
        if( (max1 < CGAL::abs(ps0)) )
        {
            max1 = CGAL::abs(ps0);
        }
        if( (max1 < CGAL::abs(pt0)) )
        {
            max1 = CGAL::abs(pt0);
        }
        if( (max1 < CGAL::abs(pu0)) )
        {
            max1 = CGAL::abs(pu0);
        }
        double max2 = CGAL::abs(pq1);
        if( (max2 < CGAL::abs(pr1)) )
        {
            max2 = CGAL::abs(pr1);
        }
        if( (max2 < CGAL::abs(ps1)) )
        {
            max2 = CGAL::abs(ps1);
        }
        if( (max2 < CGAL::abs(pt1)) )
        {
            max2 = CGAL::abs(pt1);
        }
        if( (max2 < CGAL::abs(pu1)) )
        {
            max2 = CGAL::abs(pu1);
        }
        double max3 = CGAL::abs(pq2);
        if( (max3 < CGAL::abs(pr2)) )
        {
            max3 = CGAL::abs(pr2);
        }
        if( (max3 < CGAL::abs(ps2)) )
        {
            max3 = CGAL::abs(ps2);
        }
        if( (max3 < CGAL::abs(pt2)) )
        {
            max3 = CGAL::abs(pt2);
        }
        if( (max3 < CGAL::abs(pu2)) )
        {
            max3 = CGAL::abs(pu2);
        }
        double max4 = CGAL::abs(pq3);
        if( (max4 < CGAL::abs(pr3)) )
        {
            max4 = CGAL::abs(pr3);
        }
        if( (max4 < CGAL::abs(ps3)) )
        {
            max4 = CGAL::abs(ps3);
        }
        if( (max4 < CGAL::abs(pt3)) )
        {
            max4 = CGAL::abs(pt3);
        }
        if( (max4 < CGAL::abs(pu3)) )
        {
            max4 = CGAL::abs(pu3);
        }
        double max5 = CGAL::abs(pq4);
        if( (max5 < CGAL::abs(pr4)) )
        {
            max5 = CGAL::abs(pr4);
        }
        if( (max5 < CGAL::abs(ps4)) )
        {
            max5 = CGAL::abs(ps4);
        }
        if( (max5 < CGAL::abs(pt4)) )
        {
            max5 = CGAL::abs(pt4);
        }
        if( (max5 < CGAL::abs(pu4)) )
        {
            max5 = CGAL::abs(pu4);
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
        if( (lower_bound_1 < 9.99657131447050328602e-60) )
        {
            return Base::operator()(p, q, r, s, t, u);
        }
        else
        {
            if( (upper_bound_1 > 3.21387608851797912384e+60) )
            {
                return Base::operator()(p, q, r, s, t, u);
            }
            eps = (2.22889232457534740153e-13 * ((((max1 * max2) * max3) * max4) * max5));
            if( (det > eps) )
            {
              if(should_be != POSITIVE){
                std::cout << "result should not be POSITIVE" << std::endl;
              }
                return POSITIVE;
            }

            if( (det < -eps) )
            {
              if(should_be != NEGATIVE){
                std::cout << "result should not be NEGATIVE" << std::endl;
              }
                return NEGATIVE;
            }
        }
      }
      return Base::operator()(p, q, r, s, t, u);
  }


};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_ORIENTATION_5_H
