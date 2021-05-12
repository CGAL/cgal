// Copyright (c) 2017  GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
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
    typedef typename Kernel_traits<A1>::type EK;
    typedef typename Type_mapper<A1,EK,Exact_predicates_inexact_constructions_kernel>::type T1;
    std::pair<T1,bool> aa1 = convert(a1.approx());
    if(! aa1.second){
      return fp(a1);
    }

    return epicp(aa1.first);
  }

  template <typename A1>
  result_type operator()(const A1& a1, const Null_vector& v) const
  {
    CGAL::Epic_converter<AK> convert;
    typedef typename Kernel_traits<A1>::type EK;
    typedef typename Type_mapper<A1,EK,Exact_predicates_inexact_constructions_kernel>::type T1;
    std::pair<T1,bool> aa1 = convert(a1.approx());
    if(! aa1.second){
      return fp(a1, v);
    }

    return epicp(aa1.first, v);
  }

  template <typename A1, typename A2>
  result_type operator()(const A1& a1, const A2& a2) const
  {
    CGAL::Epic_converter<AK> convert;
    typedef typename Kernel_traits<A1>::type EK;
    typedef typename Type_mapper<A1,EK, Exact_predicates_inexact_constructions_kernel>::type T1;
    std::pair<T1,bool> aa1 = convert(approx(a1));
    if(! aa1.second){
      return fp(a1, a2);
    }
    typedef typename Type_mapper<A2,EK, Exact_predicates_inexact_constructions_kernel>::type T2;
    std::pair<T2,bool> aa2 = convert(approx(a2));
    if(! aa2.second){
      return fp(a1, a2);
    }
    return epicp(aa1.first, aa2.first);
  }

    // We need these two specializations as in general we determine
    // the kernel for the template argument A1, and this does not work for Bbox_2 and Bbox_3
  template <typename A2>
  result_type operator()(const Bbox_2& bb, const A2& a2) const
  {
    CGAL::Epic_converter<AK> convert;
    typedef typename Kernel_traits<A2>::type EK;
    typedef typename Type_mapper<A2,EK, Exact_predicates_inexact_constructions_kernel>::type T2;
    std::pair<T2,bool> aa2 = convert(approx(a2));
    if(! aa2.second){
      return fp(bb, a2);
    }
    return epicp(bb, aa2.first);
  }

  template <typename A2>
  result_type operator()(const Bbox_3& bb, const A2& a2) const
  {
    CGAL::Epic_converter<AK> convert;
    typedef typename Kernel_traits<A2>::type EK;
    typedef typename Type_mapper<A2,EK, Exact_predicates_inexact_constructions_kernel>::type T2;
    std::pair<T2,bool> aa2 = convert(approx(a2));
    if(! aa2.second){
      return fp(bb, a2);
    }
    return epicp(bb, aa2.first);
  }

  template <typename A1, typename A2, typename A3>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3) const
  {
    CGAL::Epic_converter<AK> convert;
    typedef typename Kernel_traits<A1>::type EK;
    typedef typename Type_mapper<A1,EK, Exact_predicates_inexact_constructions_kernel>::type T1;
    std::pair<T1,bool> aa1 = convert(a1.approx());
    if(! aa1.second){
      return fp(a1, a2, a3);
    }
    typedef typename Type_mapper<A2, EK, Exact_predicates_inexact_constructions_kernel>::type T2;
    std::pair<T2,bool> aa2 = convert(a2.approx());
    if(! aa2.second){
      return fp(a1, a2, a3);
    }
    typedef typename Type_mapper<A3,EK, Exact_predicates_inexact_constructions_kernel>::type T3;
    std::pair<T3,bool> aa3 = convert(a3.approx());
    if(! aa3.second){
      return fp(a1, a2, a3);
    }
    return epicp(aa1.first, aa2.first, aa3.first);
  }


  template <typename A1, typename A2, typename A3, typename A4>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
  {
    CGAL::Epic_converter<AK> convert;
    typedef typename Kernel_traits<A1>::type EK;
    typedef typename Type_mapper<A1,EK,Exact_predicates_inexact_constructions_kernel>::type T1;
    std::pair<T1,bool> aa1 = convert(a1.approx());
    if(! aa1.second){
      return fp(a1, a2, a3, a4);
    }
    typedef typename Type_mapper<A2,EK,Exact_predicates_inexact_constructions_kernel>::type T2;
    std::pair<T2,bool> aa2 = convert(a2.approx());
    if(! aa2.second){
      return fp(a1, a2, a3, a4);
    }
    typedef typename Type_mapper<A3,EK,Exact_predicates_inexact_constructions_kernel>::type T3;
    std::pair<T3,bool> aa3 = convert(a3.approx());
    if(! aa3.second){
      return fp(a1, a2, a3, a4);
    }
    typedef typename Type_mapper<A4,EK,Exact_predicates_inexact_constructions_kernel>::type T4;
    std::pair<T4,bool> aa4 = convert(a4.approx());
    if(! aa4.second){
      return fp(a1, a2, a3, a4);
    }
    return epicp(aa1.first, aa2.first, aa3.first, aa4.first);
  }

  template <typename A1, typename A2, typename A3, typename A4, typename A5>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5) const
  {
    CGAL::Epic_converter<AK> convert;
    typedef typename Kernel_traits<A1>::type EK;
    typedef typename Type_mapper<A1,EK,Exact_predicates_inexact_constructions_kernel>::type T1;
    std::pair<T1,bool> aa1 = convert(a1.approx());
    if(! aa1.second){
      return fp(a1, a2, a3, a4, a5);
    }
    typedef typename Type_mapper<A2,EK,Exact_predicates_inexact_constructions_kernel>::type T2;
    std::pair<T2,bool> aa2 = convert(a2.approx());
    if(! aa2.second){
      return fp(a1, a2, a3, a4, a5);
    }
    typedef typename Type_mapper<A3,EK,Exact_predicates_inexact_constructions_kernel>::type T3;
    std::pair<T3,bool> aa3 = convert(a3.approx());
    if(! aa3.second){
      return fp(a1, a2, a3, a4, a5);
    }
    typedef typename Type_mapper<A4,EK,Exact_predicates_inexact_constructions_kernel>::type T4;
    std::pair<T4,bool> aa4 = convert(a4.approx());
    if(! aa4.second){
      return fp(a1, a2, a3, a4, a5);
    }
    typedef typename Type_mapper<A5,EK,Exact_predicates_inexact_constructions_kernel>::type T5;
    std::pair<T5,bool> aa5 = convert(a5.approx());
    if(! aa5.second){
      return fp(a1, a2, a3, a4, a5);
    }
    return epicp(aa1.first, aa2.first, aa3.first, aa4.first, aa5.first);
  }

  template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6) const
  {
    CGAL::Epic_converter<AK> convert;
    typedef typename Kernel_traits<A1>::type EK;
    typedef typename Type_mapper<A1,EK,Exact_predicates_inexact_constructions_kernel>::type T1;
    std::pair<T1,bool> aa1 = convert(a1.approx());
    if(! aa1.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    typedef typename Type_mapper<A2,EK,Exact_predicates_inexact_constructions_kernel>::type T2;
    std::pair<T2,bool> aa2 = convert(a2.approx());
    if(! aa2.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    typedef typename Type_mapper<A3,EK,Exact_predicates_inexact_constructions_kernel>::type T3;
    std::pair<T3,bool> aa3 = convert(a3.approx());
    if(! aa3.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    typedef typename Type_mapper<A4,EK,Exact_predicates_inexact_constructions_kernel>::type T4;
    std::pair<T4,bool> aa4 = convert(a4.approx());
    if(! aa4.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    typedef typename Type_mapper<A5,EK,Exact_predicates_inexact_constructions_kernel>::type T5;
    std::pair<T5,bool> aa5 = convert(a5.approx());
    if(! aa5.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    typedef typename Type_mapper<A6,EK,Exact_predicates_inexact_constructions_kernel>::type T6;
    std::pair<T6,bool> aa6 = convert(a6.approx());
    if(! aa6.second){
      return fp(a1, a2, a3, a4, a5, a6);
    }
    return epicp(aa1.first, aa2.first, aa3.first, aa4.first, aa5.first, aa6.first);
  }

  template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A6& a7) const
  {
    CGAL::Epic_converter<AK> convert;
    typedef typename Kernel_traits<A1>::type EK;
    typedef typename Type_mapper<A1,EK,Exact_predicates_inexact_constructions_kernel>::type T1;
    std::pair<T1,bool> aa1 = convert(a1.approx());
    if(! aa1.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    typedef typename Type_mapper<A2,EK,Exact_predicates_inexact_constructions_kernel>::type T2;
    std::pair<T2,bool> aa2 = convert(a2.approx());
    if(! aa2.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    typedef typename Type_mapper<A3,EK,Exact_predicates_inexact_constructions_kernel>::type T3;
    std::pair<T3,bool> aa3 = convert(a3.approx());
    if(! aa3.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    typedef typename Type_mapper<A4,EK,Exact_predicates_inexact_constructions_kernel>::type T4;
    std::pair<T4,bool> aa4 = convert(a4.approx());
    if(! aa4.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    typedef typename Type_mapper<A5,EK,Exact_predicates_inexact_constructions_kernel>::type T5;
    std::pair<T5,bool> aa5 = convert(a5.approx());
    if(! aa5.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    typedef typename Type_mapper<A6,EK,Exact_predicates_inexact_constructions_kernel>::type T6;
    std::pair<T6,bool> aa6 = convert(a6.approx());
    if(! aa6.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    typedef typename Type_mapper<A7,EK,Exact_predicates_inexact_constructions_kernel>::type T7;
    std::pair<T7,bool> aa7 = convert(a7.approx());
    if(! aa7.second){
      return fp(a1, a2, a3, a4, a5, a6, a7);
    }
    return epicp(aa1.first, aa2.first, aa3.first, aa4.first, aa5.first, aa6.first, aa7.first);
  }


  template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7, typename A8>
  result_type operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7, const A8& a8) const
  {
    CGAL::Epic_converter<AK> convert;
    typedef typename Kernel_traits<A1>::type EK;
    typedef typename Type_mapper<A1,EK,Exact_predicates_inexact_constructions_kernel>::type T1;
    std::pair<T1,bool> aa1 = convert(a1.approx());
    if(! aa1.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    typedef typename Type_mapper<A2,EK,Exact_predicates_inexact_constructions_kernel>::type T2;
    std::pair<T2,bool> aa2 = convert(a2.approx());
    if(! aa2.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    typedef typename Type_mapper<A3,EK,Exact_predicates_inexact_constructions_kernel>::type T3;
    std::pair<T3,bool> aa3 = convert(a3.approx());
    if(! aa3.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    typedef typename Type_mapper<A4,EK,Exact_predicates_inexact_constructions_kernel>::type T4;
    std::pair<T4,bool> aa4 = convert(a4.approx());
    if(! aa4.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    typedef typename Type_mapper<A5,EK,Exact_predicates_inexact_constructions_kernel>::type T5;
    std::pair<T5,bool> aa5 = convert(a5.approx());
    if(! aa5.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    typedef typename Type_mapper<A6,EK,Exact_predicates_inexact_constructions_kernel>::type T6;
    std::pair<T6,bool> aa6 = convert(a6.approx());
    if(! aa6.second){
      return fp(a1, a2, a3, a5, a5, a6, a7, a8);
    }
    typedef typename Type_mapper<A7,EK,Exact_predicates_inexact_constructions_kernel>::type T7;
    std::pair<T7,bool> aa7 = convert(a7.approx());
    if(! aa7.second){
      return fp(a1, a2, a3, a5, a5, a6, a7, a8);
    }
    typedef typename Type_mapper<A8,EK,Exact_predicates_inexact_constructions_kernel>::type T8;
    std::pair<T8,bool> aa8 = convert(a8.approx());
    if(! aa8.second){
      return fp(a1, a2, a3, a4, a5, a6, a7, a8);
    }
    return epicp(aa1.first, aa2.first, aa3.first, aa4.first, aa5.first, aa6.first, aa7.first, aa8.first);
  }
};

} // CGAL

#endif // CGAL_STATIC_FILTERED_PREDICATE_H
