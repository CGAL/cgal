// Copyright (c) 2006-2008  Max-Planck-Institute Saarbruecken (Germany)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Michael Hemmer <hemmer@mpi-sb.mpg.de>

#ifndef CGAL_UTILS_CLASSES_H
#define CGAL_UTILS_CLASSES_H
#include <CGAL/basic.h>

namespace CGAL {

template < class A, class B = A >
struct Equal_to : public std::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x == y; }
};

template < class A, class B = A >
struct Not_equal_to : public std::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x != y; }
};

template < class A, class B = A >
struct Greater : public std::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x > y; }
};

template < class A, class B = A >
struct Less : public std::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x < y; }
};

template < class A, class B = A >
struct Greater_equal : public std::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x >= y; }
};

template < class A, class B = A >
struct Less_equal : public std::binary_function< A, B, bool > {
  bool operator()( const A& x, const B& y) const
  { return x <= y; }
};

template < class NT, class Less = std::less< NT > >
struct Min :public std::binary_function< NT, NT, NT > {
 Min() {}
 Min(const Less& c_) : c(c_) {}
 NT operator()( const NT& x, const NT& y) const
    { return (std::min) BOOST_PREVENT_MACRO_SUBSTITUTION ( x, y, c); }
protected:
 Less c;
};

template < class NT, class Less = std::less< NT > >
struct Max :public std::binary_function< NT, NT, NT > {
 Max() {}
 Max(const Less& c_) : c(c_) {}
 NT operator()( const NT& x, const NT& y) const
    { return (std::max) BOOST_PREVENT_MACRO_SUBSTITUTION ( x, y, c); }
protected:
 Less c;
};

template< class T >
class Is_valid
  : public std::unary_function< T, bool > {
  public:
    bool operator()( const T& ) const {
      return true;
    };
};

} //namespace CGAL

#endif // CGAL_UTILS_CLASSES_H
