// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_ARRAY_H
#define CGAL_ARRAY_H

#include <CGAL/config.h>
#ifndef CGAL_CFG_NO_CPP0X_ARRAY
#  include <array>
#else
#  include <boost/array.hpp>
#endif

namespace CGAL {

// The make_array() function simply constructs an std::array.
// It is needed for cases where a std::array is used as a class data
// member and you want to initialize it in the member initializers list.
// It is also optimized: no spurious copies of the objects are made,
// provided the compiler does the NRVO.  So this is better than
// default construction followed by assignment.

// I proposed it for Boost.Array, but it has not been integrated so far.
// See the thread starting at
// https://lists.boost.org/Archives/boost/2006/08/109003.php
//
// C++0x has it under discussion here :
// http://www.open-std.org/jtc1/sc22/wg21/docs/lwg-active.html#851

// Hopefully C++0x will fix this properly with initializer_lists.
// So, it's temporary, therefore I do not document it and keep it internal.

// NOTE : The above is actually untrue !  It is possible to do :
//     struct S2 {
//       typedef boost::array<M,2> Ar;
//       Ar m;
//       S2 (const M&a) : m ((Ar) { { a, a } }) {}
//     };
// without spurious copies...  Except VC++ does not eat it :-(

// It's also untrue that this is not documented...  It is !

template< typename T, typename... Args >
inline
std::array< T, 1 + sizeof...(Args) >
make_array(const T & t, const Args & ... args)
{
  std::array< T, 1 + sizeof...(Args) > a = { { t, static_cast<T>(args)... } };
  return a;
}


// Functor version
struct Construct_array
{
  template <typename T, typename... Args>
  std::array<T, 1 + sizeof...(Args)> operator()(const T& t, const Args& ... args)
  {
    return make_array (t, args...);
  }
};

} //namespace CGAL

#endif // CGAL_ARRAY_H
