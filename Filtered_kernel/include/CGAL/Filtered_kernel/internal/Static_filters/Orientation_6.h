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

      double a01, a02, a03, a04, a05, a06, a11, a12, a13, a14, a15, a16,
             a21, a22, a23, a24, a25, a26, a31, a32, a33, a34, a35,
             a36, a41, a42, a43, a44, a45, a46,
             a51, a52, a53, a54, a55, a56,
             a61, a62, a63, a64, a65, a66;

      if (fit_in_double(p.c0(), a01) && fit_in_double(p.c1(), a02) &&
          fit_in_double(p.c2(), a03) && fit_in_double(p.c3(), a04) &&
          fit_in_double(p.c4(), a05) && fit_in_double(p.c5(), a06) &&

          fit_in_double(q.c0(), a11) && fit_in_double(q.c1(), a12) &&
          fit_in_double(q.c2(), a13) && fit_in_double(q.c3(), a14) &&
          fit_in_double(q.c4(), a15) && fit_in_double(q.c5(), a16) &&

          fit_in_double(r.c0(), a21) && fit_in_double(r.c1(), a22) &&
          fit_in_double(r.c2(), a23) && fit_in_double(r.c3(), a24) &&
          fit_in_double(r.c4(), a25) && fit_in_double(r.c5(), a26) &&

          fit_in_double(s.c0(), a31) && fit_in_double(s.c1(), a32) &&
          fit_in_double(s.c2(), a33) && fit_in_double(s.c3(), a34) &&
          fit_in_double(s.c4(), a35) && fit_in_double(s.c5(), a36) &&

          fit_in_double(t.c0(), a41) && fit_in_double(t.c1(), a42) &&
          fit_in_double(t.c2(), a43) && fit_in_double(t.c3(), a44) &&
          fit_in_double(t.c4(), a45) && fit_in_double(t.c5(), a46) &&

          fit_in_double(u.c0(), a51) && fit_in_double(u.c1(), a52) &&
          fit_in_double(u.c2(), a53) && fit_in_double(u.c3(), a54) &&
          fit_in_double(u.c4(), a55) && fit_in_double(u.c5(), a56) &&

          fit_in_double(v.c0(), a61) && fit_in_double(v.c1(), a62) &&
          fit_in_double(v.c2(), a63) && fit_in_double(v.c3(), a64) &&
          fit_in_double(v.c4(), a65) && fit_in_double(v.c5(), a66))
      {
        CGAL_assertion_code(Orientation should_be = Base::operator()(p, q, r, s, t, u, v));
        double m01;
    m01 = (a11 - a01);
    double m02;
    m02 = (a21 - a01);
    double m03;
    m03 = (a31 - a01);
    double m04;
    m04 = (a41 - a01);
    double m05;
    m05 = (a51 - a01);
    double m06;
    m06 = (a61 - a01);
    double m12;
    m12 = (a21 - a11);
    double m13;
    m13 = (a31 - a11);
    double m14;
    m14 = (a41 - a11);
    double m15;
    m15 = (a51 - a11);
    double m16;
    m16 = (a61 - a11);
    double m23;
    m23 = (a31 - a21);
    double m24;
    m24 = (a41 - a21);
    double m25;
    m25 = (a51 - a21);
    double m26;
    m26 = (a61 - a21);
    double m34;
    m34 = (a41 - a31);
    double m35;
    m35 = (a51 - a31);
    double m36;
    m36 = (a61 - a31);
    double m45;
    m45 = (a51 - a41);
    double m46;
    m46 = (a61 - a41);
    double m56;
    m56 = (a61 - a51);
    double m012;
    m012 = (((m01 * a22) - (m02 * a12)) + (m12 * a02));
    double m013;
    m013 = (((m01 * a32) - (m03 * a12)) + (m13 * a02));
    double m014;
    m014 = (((m01 * a42) - (m04 * a12)) + (m14 * a02));
    double m015;
    m015 = (((m01 * a52) - (m05 * a12)) + (m15 * a02));
    double m016;
    m016 = (((m01 * a62) - (m06 * a12)) + (m16 * a02));
    double m023;
    m023 = (((m02 * a32) - (m03 * a22)) + (m23 * a02));
    double m024;
    m024 = (((m02 * a42) - (m04 * a22)) + (m24 * a02));
    double m025;
    m025 = (((m02 * a52) - (m05 * a22)) + (m25 * a02));
    double m026;
    m026 = (((m02 * a62) - (m06 * a22)) + (m26 * a02));
    double m034;
    m034 = (((m03 * a42) - (m04 * a32)) + (m34 * a02));
    double m035;
    m035 = (((m03 * a52) - (m05 * a32)) + (m35 * a02));
    double m036;
    m036 = (((m03 * a62) - (m06 * a32)) + (m36 * a02));
    double m045;
    m045 = (((m04 * a52) - (m05 * a42)) + (m45 * a02));
    double m046;
    m046 = (((m04 * a62) - (m06 * a42)) + (m46 * a02));
    double m056;
    m056 = (((m05 * a62) - (m06 * a52)) + (m56 * a02));
    double m123;
    m123 = (((m12 * a32) - (m13 * a22)) + (m23 * a12));
    double m124;
    m124 = (((m12 * a42) - (m14 * a22)) + (m24 * a12));
    double m125;
    m125 = (((m12 * a52) - (m15 * a22)) + (m25 * a12));
    double m126;
    m126 = (((m12 * a62) - (m16 * a22)) + (m26 * a12));
    double m134;
    m134 = (((m13 * a42) - (m14 * a32)) + (m34 * a12));
    double m135;
    m135 = (((m13 * a52) - (m15 * a32)) + (m35 * a12));
    double m136;
    m136 = (((m13 * a62) - (m16 * a32)) + (m36 * a12));
    double m145;
    m145 = (((m14 * a52) - (m15 * a42)) + (m45 * a12));
    double m146;
    m146 = (((m14 * a62) - (m16 * a42)) + (m46 * a12));
    double m156;
    m156 = (((m15 * a62) - (m16 * a52)) + (m56 * a12));
    double m234;
    m234 = (((m23 * a42) - (m24 * a32)) + (m34 * a22));
    double m235;
    m235 = (((m23 * a52) - (m25 * a32)) + (m35 * a22));
    double m236;
    m236 = (((m23 * a62) - (m26 * a32)) + (m36 * a22));
    double m245;
    m245 = (((m24 * a52) - (m25 * a42)) + (m45 * a22));
    double m246;
    m246 = (((m24 * a62) - (m26 * a42)) + (m46 * a22));
    double m256;
    m256 = (((m25 * a62) - (m26 * a52)) + (m56 * a22));
    double m345;
    m345 = (((m34 * a52) - (m35 * a42)) + (m45 * a32));
    double m346;
    m346 = (((m34 * a62) - (m36 * a42)) + (m46 * a32));
    double m356;
    m356 = (((m35 * a62) - (m36 * a52)) + (m56 * a32));
    double m456;
    m456 = (((m45 * a62) - (m46 * a52)) + (m56 * a42));
    double m0123;
    m0123 = ((((m012 * a33) - (m013 * a23)) + (m023 * a13)) - (m123 * a03));
    double m0124;
    m0124 = ((((m012 * a43) - (m014 * a23)) + (m024 * a13)) - (m124 * a03));
    double m0125;
    m0125 = ((((m012 * a53) - (m015 * a23)) + (m025 * a13)) - (m125 * a03));
    double m0126;
    m0126 = ((((m012 * a63) - (m016 * a23)) + (m026 * a13)) - (m126 * a03));
    double m0134;
    m0134 = ((((m013 * a43) - (m014 * a33)) + (m034 * a13)) - (m134 * a03));
    double m0135;
    m0135 = ((((m013 * a53) - (m015 * a33)) + (m035 * a13)) - (m135 * a03));
    double m0136;
    m0136 = ((((m013 * a63) - (m016 * a33)) + (m036 * a13)) - (m136 * a03));
    double m0145;
    m0145 = ((((m014 * a53) - (m015 * a43)) + (m045 * a13)) - (m145 * a03));
    double m0146;
    m0146 = ((((m014 * a63) - (m016 * a43)) + (m046 * a13)) - (m146 * a03));
    double m0156;
    m0156 = ((((m015 * a63) - (m016 * a53)) + (m056 * a13)) - (m156 * a03));
    double m0234;
    m0234 = ((((m023 * a43) - (m024 * a33)) + (m034 * a23)) - (m234 * a03));
    double m0235;
    m0235 = ((((m023 * a53) - (m025 * a33)) + (m035 * a23)) - (m235 * a03));
    double m0236;
    m0236 = ((((m023 * a63) - (m026 * a33)) + (m036 * a23)) - (m236 * a03));
    double m0245;
    m0245 = ((((m024 * a53) - (m025 * a43)) + (m045 * a23)) - (m245 * a03));
    double m0246;
    m0246 = ((((m024 * a63) - (m026 * a43)) + (m046 * a23)) - (m246 * a03));
    double m0256;
    m0256 = ((((m025 * a63) - (m026 * a53)) + (m056 * a23)) - (m256 * a03));
    double m0345;
    m0345 = ((((m034 * a53) - (m035 * a43)) + (m045 * a33)) - (m345 * a03));
    double m0346;
    m0346 = ((((m034 * a63) - (m036 * a43)) + (m046 * a33)) - (m346 * a03));
    double m0356;
    m0356 = ((((m035 * a63) - (m036 * a53)) + (m056 * a33)) - (m356 * a03));
    double m0456;
    m0456 = ((((m045 * a63) - (m046 * a53)) + (m056 * a43)) - (m456 * a03));
    double m1234;
    m1234 = ((((m123 * a43) - (m124 * a33)) + (m134 * a23)) - (m234 * a13));
    double m1235;
    m1235 = ((((m123 * a53) - (m125 * a33)) + (m135 * a23)) - (m235 * a13));
    double m1236;
    m1236 = ((((m123 * a63) - (m126 * a33)) + (m136 * a23)) - (m236 * a13));
    double m1245;
    m1245 = ((((m124 * a53) - (m125 * a43)) + (m145 * a23)) - (m245 * a13));
    double m1246;
    m1246 = ((((m124 * a63) - (m126 * a43)) + (m146 * a23)) - (m246 * a13));
    double m1256;
    m1256 = ((((m125 * a63) - (m126 * a53)) + (m156 * a23)) - (m256 * a13));
    double m1345;
    m1345 = ((((m134 * a53) - (m135 * a43)) + (m145 * a33)) - (m345 * a13));
    double m1346;
    m1346 = ((((m134 * a63) - (m136 * a43)) + (m146 * a33)) - (m346 * a13));
    double m1356;
    m1356 = ((((m135 * a63) - (m136 * a53)) + (m156 * a33)) - (m356 * a13));
    double m1456;
    m1456 = ((((m145 * a63) - (m146 * a53)) + (m156 * a43)) - (m456 * a13));
    double m2345;
    m2345 = ((((m234 * a53) - (m235 * a43)) + (m245 * a33)) - (m345 * a23));
    double m2346;
    m2346 = ((((m234 * a63) - (m236 * a43)) + (m246 * a33)) - (m346 * a23));
    double m2356;
    m2356 = ((((m235 * a63) - (m236 * a53)) + (m256 * a33)) - (m356 * a23));
    double m2456;
    m2456 = ((((m245 * a63) - (m246 * a53)) + (m256 * a43)) - (m456 * a23));
    double m3456;
    m3456 = ((((m345 * a63) - (m346 * a53)) + (m356 * a43)) - (m456 * a33));
    double m01234;
    m01234 = (((((m0123 * a44) - (m0124 * a34)) + (m0134 * a24)) - (m0234 * a14)) + (m1234 * a04));
    double m01235;
    m01235 = (((((m0123 * a54) - (m0125 * a34)) + (m0135 * a24)) - (m0235 * a14)) + (m1235 * a04));
    double m01236;
    m01236 = (((((m0123 * a64) - (m0126 * a34)) + (m0136 * a24)) - (m0236 * a14)) + (m1236 * a04));
    double m01245;
    m01245 = (((((m0124 * a54) - (m0125 * a44)) + (m0145 * a24)) - (m0245 * a14)) + (m1245 * a04));
    double m01246;
    m01246 = (((((m0124 * a64) - (m0126 * a44)) + (m0146 * a24)) - (m0246 * a14)) + (m1246 * a04));
    double m01256;
    m01256 = (((((m0125 * a64) - (m0126 * a54)) + (m0156 * a24)) - (m0256 * a14)) + (m1256 * a04));
    double m01345;
    m01345 = (((((m0134 * a54) - (m0135 * a44)) + (m0145 * a34)) - (m0345 * a14)) + (m1345 * a04));
    double m01346;
    m01346 = (((((m0134 * a64) - (m0136 * a44)) + (m0146 * a34)) - (m0346 * a14)) + (m1346 * a04));
    double m01356;
    m01356 = (((((m0135 * a64) - (m0136 * a54)) + (m0156 * a34)) - (m0356 * a14)) + (m1356 * a04));
    double m01456;
    m01456 = (((((m0145 * a64) - (m0146 * a54)) + (m0156 * a44)) - (m0456 * a14)) + (m1456 * a04));
    double m02345;
    m02345 = (((((m0234 * a54) - (m0235 * a44)) + (m0245 * a34)) - (m0345 * a24)) + (m2345 * a04));
    double m02346;
    m02346 = (((((m0234 * a64) - (m0236 * a44)) + (m0246 * a34)) - (m0346 * a24)) + (m2346 * a04));
    double m02356;
    m02356 = (((((m0235 * a64) - (m0236 * a54)) + (m0256 * a34)) - (m0356 * a24)) + (m2356 * a04));
    double m02456;
    m02456 = (((((m0245 * a64) - (m0246 * a54)) + (m0256 * a44)) - (m0456 * a24)) + (m2456 * a04));
    double m03456;
    m03456 = (((((m0345 * a64) - (m0346 * a54)) + (m0356 * a44)) - (m0456 * a34)) + (m3456 * a04));
    double m12345;
    m12345 = (((((m1234 * a54) - (m1235 * a44)) + (m1245 * a34)) - (m1345 * a24)) + (m2345 * a14));
    double m12346;
    m12346 = (((((m1234 * a64) - (m1236 * a44)) + (m1246 * a34)) - (m1346 * a24)) + (m2346 * a14));
    double m12356;
    m12356 = (((((m1235 * a64) - (m1236 * a54)) + (m1256 * a34)) - (m1356 * a24)) + (m2356 * a14));
    double m12456;
    m12456 = (((((m1245 * a64) - (m1246 * a54)) + (m1256 * a44)) - (m1456 * a24)) + (m2456 * a14));
    double m13456;
    m13456 = (((((m1345 * a64) - (m1346 * a54)) + (m1356 * a44)) - (m1456 * a34)) + (m3456 * a14));
    double m23456;
    m23456 = (((((m2345 * a64) - (m2346 * a54)) + (m2356 * a44)) - (m2456 * a34)) + (m3456 * a24));
    double m012345;
    m012345 = ((((((m01234 * a55) - (m01235 * a45)) + (m01245 * a35)) - (m01345 * a25)) + (m02345 * a15)) - (m12345 * a05));
    double m012346;
    m012346 = ((((((m01234 * a65) - (m01236 * a45)) + (m01246 * a35)) - (m01346 * a25)) + (m02346 * a15)) - (m12346 * a05));
    double m012356;
    m012356 = ((((((m01235 * a65) - (m01236 * a55)) + (m01256 * a35)) - (m01356 * a25)) + (m02356 * a15)) - (m12356 * a05));
    double m012456;
    m012456 = ((((((m01245 * a65) - (m01246 * a55)) + (m01256 * a45)) - (m01456 * a25)) + (m02456 * a15)) - (m12456 * a05));
    double m013456;
    m013456 = ((((((m01345 * a65) - (m01346 * a55)) + (m01356 * a45)) - (m01456 * a35)) + (m03456 * a15)) - (m13456 * a05));
    double m023456;
    m023456 = ((((((m02345 * a65) - (m02346 * a55)) + (m02356 * a45)) - (m02456 * a35)) + (m03456 * a25)) - (m23456 * a05));
    double m123456;
    m123456 = ((((((m12345 * a65) - (m12346 * a55)) + (m12356 * a45)) - (m12456 * a35)) + (m13456 * a25)) - (m23456 * a15));
    double m0123456;
    m0123456 = (((((((m012345 * a66) - (m012346 * a56)) + (m012356 * a46)) - (m012456 * a36)) + (m013456 * a26)) - (m023456 * a16)) + (m123456 * a06));
    int int_tmp_result;
    double eps;
    double max1 = CGAL::abs(a02);
    double tmp = CGAL::abs(a12);
    if( max1 < tmp )
    {
        max1 = tmp;
    }
    tmp = CGAL::abs(a22);
    if( max1 < tmp )
    {
        max1 = tmp;
    }
    tmp = CGAL::abs(a32);
    if( max1 < tmp)
    {
        max1 = tmp;
    }
    tmp = CGAL::abs(a42);
    if( max1 < tmp)
    {
        max1 = tmp;
    }
    tmp = CGAL::abs(a52);
    if(max1 < tmp )
    {
        max1 = tmp;
    }
    tmp = CGAL::abs(a62);
    if( (max1 < tmp) )
    {
        max1 = tmp;
    }
    double max2 = CGAL::abs(a03);
    tmp = CGAL::abs(a13);
    if( (max2 < tmp) )
    {
        max2 = tmp;
    }
    tmp = CGAL::abs(a23);
    if( (max2 < tmp) )
    {
        max2 = tmp;
    }
    tmp = CGAL::abs(a33);
    if( (max2 < tmp) )
    {
        max2 = tmp;
    }
    tmp = CGAL::abs(a43);
    if( (max2 < tmp) )
    {
        max2 = tmp;
    }
    tmp = CGAL::abs(a53);
    if( (max2 < tmp) )
    {
        max2 = tmp;
    }
    tmp = CGAL::abs(a63);
    if( (max2 < tmp) )
    {
        max2 = tmp;
    }
    double max3 = CGAL::abs(a04);
    tmp = CGAL::abs(a14);
    if( (max3 < tmp) )
    {
        max3 = tmp;
    }
    tmp = CGAL::abs(a24);
    if( (max3 < tmp) )
    {
        max3 = tmp;
    }
    tmp = CGAL::abs(a34);
    if( (max3 < tmp) )
    {
        max3 = tmp;
    }
    tmp = CGAL::abs(a44);
    if( (max3 < tmp) )
    {
        max3 = tmp;
    }
    tmp = CGAL::abs(a54);
    if( (max3 < tmp) )
    {
        max3 = tmp;
    }
    tmp = CGAL::abs(a64);
    if( (max3 < tmp) )
    {
        max3 = tmp;
    }
    double max4 = CGAL::abs(a05);
    tmp = CGAL::abs(a15);
    if( (max4 < tmp) )
    {
        max4 = tmp;
    }
    tmp = CGAL::abs(a25);
    if( (max4 < tmp) )
    {
        max4 = tmp;
    }
    tmp = CGAL::abs(a35);
    if( (max4 < tmp) )
    {
        max4 = tmp;
    }
    tmp = CGAL::abs(a45);
    if( (max4 < tmp) )
    {
        max4 = tmp;
    }
    tmp = CGAL::abs(a55);
    if( (max4 < tmp) )
    {
        max4 = tmp;
    }
    tmp = CGAL::abs(a65);
    if( (max4 < tmp) )
    {
        max4 = tmp;
    }
    double max5 = CGAL::abs(a06);
    tmp = CGAL::abs(a16);
    if( (max5 < tmp) )
    {
        max5 = tmp;
    }
    tmp = CGAL::abs(a26);
    if( (max5 < tmp) )
    {
        max5 = tmp;
    }
    tmp = CGAL::abs(a36);
    if( (max5 < tmp) )
    {
        max5 = tmp;
    }
    tmp = CGAL::abs(a46);
    if( (max5 < tmp) )
    {
        max5 = tmp;
    }
    tmp = CGAL::abs(a56);
    if( (max5 < tmp) )
    {
        max5 = tmp;
    }
    tmp = CGAL::abs(a66);
    if( (max5 < tmp) )
    {
        max5 = tmp;
    }

    double max6 = CGAL::abs(m01);
    tmp = CGAL::abs(m02);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m03);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m04);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m05);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m06);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m12);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m13);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m14);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m15);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m16);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m23);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m24);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m25);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m26);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m34);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m35);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m36);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m45);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m46);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }
    tmp = CGAL::abs(m56);
    if( (max6 < tmp) )
    {
        max6 = tmp;
    }

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
    if( (lower_bound_1 < 3.98278627399642140002e-50) )
    {
        return Base::operator()(p, q, r, s, t, u, v);
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355511223e+50) )
        {
            return Base::operator()(p, q, r, s, t, u, v);
        }
        eps = (5.57471180948088246730e-12 * (((((max6 * max1) * max2) * max3) * max4) * max5));
        if( (m0123456 > eps) )
        {
            CGAL_assertion(should_be == POSITIVE);
            return POSITIVE;
        }
        else
        {
            if( (m0123456 < -eps) )
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
