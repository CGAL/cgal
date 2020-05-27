// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_TYPESET_H
#define CGAL_TYPESET_H
#include <CGAL/config.h>
#include <type_traits>

// Sometimes using tuple just to list types is overkill (takes forever to
// instantiate).

namespace CGAL {
  template<class...> struct typeset;
  template<class H,class...U> struct typeset<H,U...> {
    typedef H head;
    typedef typeset<U...> tail;
    typedef typeset type;
    template<class X> using contains = typename
      std::conditional<
        std::is_same<H,X>::value,
        std::true_type,
        typename tail::template contains<X>
      >::type;
    template<class X> using add = typename
      std::conditional<
        contains<X>::value,
        typeset<H,U...>,
        typeset<H,U...,X>
      >::type;
  };
  template<> struct typeset<> {
    typedef typeset type;
    template<class X> using contains = std::false_type;
    template<class X> using add = typeset<X>;
  };
  struct typeset_all {
    typedef typeset_all type;
    template<class X> using contains = std::true_type;
    template<class X> using add = typeset_all;
  };

  template<class T1, class T2> struct typeset_union_ :
    typeset_union_<typename T1::template add<typename T2::head>::type, typename T2::tail>
  {};
  template<class T> struct typeset_union_ <T, typeset<> > : T {};
  template<class T> struct typeset_union_ <T, typeset_all > : typeset_all {};

  template<class T1, class T2>
    struct typeset_intersection_ {
      typedef typename T1::head H;
      typedef typename typeset_intersection_<typename T1::tail,T2>::type U;
      typedef typename
        std::conditional<T2::template contains<H>::value,
        typename U::template add<H>::type, U>::type type;
    };
  template<class T> struct typeset_intersection_<typeset<>, T> : typeset<> {};
  template<class T> struct typeset_intersection_<typeset_all, T> : T {};

  template<class T1, class T2>
    using typeset_union = typename typeset_union_<T1,T2>::type;
  template<class T1, class T2>
    using typeset_intersection = typename typeset_intersection_<T1,T2>::type;
}
#endif
