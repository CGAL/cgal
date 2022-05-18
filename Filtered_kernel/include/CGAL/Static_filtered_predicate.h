// Copyright (c) 2017  GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri, Laurent Rineau

#ifndef CGAL_STATIC_FILTERED_PREDICATE_H
#define CGAL_STATIC_FILTERED_PREDICATE_H

#include <CGAL/Epic_converter.h>

namespace CGAL {

template <typename AK, typename FP, typename EpicP>
class Static_filtered_predicate {
public:
  FP fp;
  EpicP epicp;
  typedef typename AK::FT IA;
  typedef typename FP::result_type result_type;

  template <typename A1>
  result_type operator()(const A1& a1) const
  {
    CGAL::Epic_converter<AK> convert;
    auto aa1 = convert(approx(a1));
    if(! aa1.second){
      return fp(a1);
    }

    return epicp(aa1.first);
  }


  template <typename A1, typename A2>
  result_type operator()(const A1& a1, const A2& a2) const
  {
    CGAL::Epic_converter<AK> convert;
    auto aa1 = convert(approx(a1));
    if(! aa1.second){
      return fp(a1, a2);
    }
    auto aa2 = convert(approx(a2));
    if(! aa2.second){
      return fp(a1, a2);
    }
    return epicp(aa1.first, aa2.first);
  }


  template <typename A1, typename A2, typename A3>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3) const
  {
    CGAL::Epic_converter<AK> convert;
    auto aa1 = convert(approx(a1));
    if(! aa1.second){
      return fp(a1, a2, a3);
    }
    auto aa2 = convert(approx(a2));
    if(! aa2.second){
      return fp(a1, a2, a3);
    }
    auto aa3 = convert(approx(a3));
    if(! aa3.second){
      return fp(a1, a2, a3);
    }
    return epicp(aa1.first, aa2.first, aa3.first);
  }


  template <typename A1, typename A2, typename A3, typename A4>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
  {
    CGAL::Epic_converter<AK> convert;
    auto aa1 = convert(approx(a1));
    if(! aa1.second){
      return fp(a1, a2, a3, a4);
    }

    auto aa2 = convert(approx(a2));
    if(! aa2.second){
      return fp(a1, a2, a3, a4);
    }

    auto aa3 = convert(approx(a3));
    if(! aa3.second){
      return fp(a1, a2, a3, a4);
    }

    auto aa4 = convert(approx(a4));
    if(! aa4.second){
      return fp(a1, a2, a3, a4);
    }
    return epicp(aa1.first, aa2.first, aa3.first, aa4.first);
  }

  template <typename A1, typename A2, typename A3, typename A4, typename A5>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5) const
  {
    CGAL::Epic_converter<AK> convert;
    auto aa1 = convert(approx(a1));
    if(! aa1.second){
      return fp(a1, a2, a3, a4, a5);
    }
    auto aa2 = convert(approx(a2));
    if(! aa2.second){
      return fp(a1, a2, a3, a4, a5);
    }
    auto aa3 = convert(approx(a3));
    if(! aa3.second){
      return fp(a1, a2, a3, a4, a5);
    }
    auto aa4 = convert(approx(a4));
    if(! aa4.second){
      return fp(a1, a2, a3, a4, a5);
    }
    auto aa5 = convert(approx(a5));
    if(! aa5.second){
      return fp(a1, a2, a3, a4, a5);
    }
    return epicp(aa1.first, aa2.first, aa3.first, aa4.first, aa5.first);
  }

  template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6) const
  {
    CGAL::Epic_converter<AK> convert;
    auto aa1 = convert(approx(a1));
    if(! aa1.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    auto aa2 = convert(approx(a2));
    if(! aa2.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    auto aa3 = convert(approx(a3));
    if(! aa3.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    auto aa4 = convert(approx(a4));
    if(! aa4.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    auto aa5 = convert(approx(a5));
    if(! aa5.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    auto aa6 = convert(approx(a6));
    if(! aa6.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    return epicp(aa1.first, aa2.first, aa3.first, aa4.first, aa5.first, aa6.first);
  }

  template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A6& a7) const
  {
    CGAL::Epic_converter<AK> convert;
    auto aa1 = convert(approx(a1));
    if(! aa1.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    auto aa2 = convert(approx(a2));
    if(! aa2.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    auto aa3 = convert(approx(a3));
    if(! aa3.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    auto aa4 = convert(approx(a4));
    if(! aa4.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    auto aa5 = convert(approx(a5));
    if(! aa5.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    auto aa6 = convert(approx(a6));
    if(! aa6.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    auto aa7 = convert(approx(a7));
    if(! aa7.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    return epicp(aa1.first, aa2.first, aa3.first, aa4.first, aa5.first, aa6.first, aa7.first);
  }


  template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7, typename A8>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7, const A8& a8) const
  {
    CGAL::Epic_converter<AK> convert;
    auto aa1 = convert(approx(a1));
    if(! aa1.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    auto aa2 = convert(approx(a2));
    if(! aa2.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    auto aa3 = convert(approx(a3));
    if(! aa3.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    auto aa4 = convert(approx(a4));
    if(! aa4.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    auto aa5 = convert(approx(a5));
    if(! aa5.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    auto aa6 = convert(approx(a6));
    if(! aa6.second){
      return fp(a1, a2, a3, a5, a5, a6, a7, a8);
    }
    auto aa7 = convert(approx(a7));
    if(! aa7.second){
      return fp(a1, a2, a3, a5, a5, a6, a7, a8);
    }
    auto aa8 = convert(approx(a8));
    if(! aa8.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    return epicp(aa1.first, aa2.first, aa3.first, aa4.first, aa5.first, aa6.first, aa7.first, aa8.first);
  }
};

} // CGAL

#endif // CGAL_STATIC_FILTERED_PREDICATE_H
