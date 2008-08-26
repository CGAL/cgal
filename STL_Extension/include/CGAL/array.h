// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/STL_Extension/include/CGAL/Fourtuple.h $
// $Id: Fourtuple.h 28567 2006-02-16 14:30:13Z lsaboret $
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_ARRAY_H
#define CGAL_ARRAY_H

#include <CGAL/config.h>
#ifndef CGAL_CFG_NO_CPP0X_ARRAY
#  include <array>
#elif !defined CGAL_CFG_NO_TR1_ARRAY
#  include <tr1/array>
#else
#  include <boost/array.hpp>
#endif

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_CPP0X_ARRAY
using std::array;
#elif !defined CGAL_CFG_NO_TR1_ARRAY
using std::tr1::array;
#else
using boost::array;
#endif

// The make_array() function simply constructs an std::array.
// It is needed for cases where a std::array is used as a class data
// member and you want to initialize it in the member initializers list.
// It is also optimized: no spurious copies of the objects are made,
// provided the compiler does the NRVO.  So this is better than
// default construction followed by assignment.

// I proposed it for Boost.Array, but it has not been integrated so far.
// See the thread starting at
// http://lists.boost.org/Archives/boost/2006/08/109003.php

// Hopefully C++0x will fix this properly with initializer_lists.
// So, it's temporary, therefore I do not document it and keep it internal.

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

template< typename T, typename... Args >
inline
array< T, 1 + sizeof...(Args) >
make_array(const T & t, const Args & ... args)
{
  array< T, 1 + sizeof...(Args) > a = { { t, args... } };
  return a;
}

#else // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

template < typename T > inline
array<T, 1>
make_array(const T& b1)
{
  array<T, 1> a = { { b1 } };
  return a;
}

template < typename T > inline
array<T, 2>
make_array(const T& b1, const T& b2)
{
  array<T, 2> a = { { b1, b2 } };
  return a;
}

template < typename T > inline
array<T, 3>
make_array(const T& b1, const T& b2, const T& b3)
{
  array<T, 3> a = { { b1, b2, b3 } };
  return a;
}

template < typename T > inline
array<T, 4>
make_array(const T& b1, const T& b2, const T& b3, const T& b4)
{
  array<T, 4> a = { { b1, b2, b3, b4 } };
  return a;
}

template < typename T > inline
array<T, 5>
make_array(const T& b1, const T& b2, const T& b3, const T& b4, const T& b5)
{
  array<T, 5> a = { { b1, b2, b3, b4, b5 } };
  return a;
}

template < typename T > inline
array<T, 6>
make_array(const T& b1, const T& b2, const T& b3, const T& b4, const T& b5,
           const T& b6)
{
  array<T, 6> a = { { b1, b2, b3, b4, b5, b6 } };
  return a;
}

#endif // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

CGAL_END_NAMESPACE

#endif // CGAL_ARRAY_H
