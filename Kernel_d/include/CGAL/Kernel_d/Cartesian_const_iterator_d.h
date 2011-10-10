// Copyright (c) 2000,2001,2008  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Michael Seel, Sylvain Pion

#ifndef CGAL_CARTESIAN_CONST_ITERATOR_D_H
#define CGAL_CARTESIAN_CONST_ITERATOR_D_H

#include <CGAL/basic.h>
#include <CGAL/Quotient.h>
#include <iterator>

namespace CGAL {

// Takes an iterator over RT, and make one over FT, by dividing
// by the last element.

template < typename RT_iterator >
class Cartesian_const_iterator_d
{
  typedef typename std::iterator_traits<RT_iterator>::value_type RT;
  typedef Cartesian_const_iterator_d      self;

public:

  typedef std::random_access_iterator_tag iterator_category;
  typedef CGAL::Quotient<RT>              value_type;
  typedef std::ptrdiff_t                  difference_type;
  typedef const value_type*               pointer;
  typedef const value_type&               reference;

  Cartesian_const_iterator_d() {}
  Cartesian_const_iterator_d(RT_iterator it, RT_iterator w)
    : _it(it), _w(w) {}

  self& operator++() { ++_it; return *this; }
  self  operator++(int) { self tmp = *this; ++_it; return tmp; }
  self& operator--() { --_it; return *this; }
  self  operator--(int) { self tmp = *this; --_it; return tmp; }

  self& operator+=(difference_type i) { _it+=i; return *this; }
  self& operator-=(difference_type i) { _it-=i; return *this; }
  self operator+(difference_type i) const
  { self tmp=*this; return tmp += i; }
  self operator-(difference_type i) const
  { self tmp=*this; return tmp -= i; }

  difference_type operator-(self x) const { return _it-x._it; }

  value_type operator*() const { return value_type(*_it,*_w); }
  value_type operator[](difference_type i) const { return *(*this + i); }

  bool operator==(const self& x) const { return _it==x._it; }
  bool operator!=(const self& x) const { return ! (*this==x); }
  bool operator<(const self& x) const { return (x - *this) > 0; }

private:
  RT_iterator _it, _w;
};

template < typename RT_iterator > inline
Cartesian_const_iterator_d<RT_iterator>
operator+(std::ptrdiff_t i, Cartesian_const_iterator_d<RT_iterator> const& it)
{
	return it+i;
}

template < typename RT_iterator > inline
Cartesian_const_iterator_d<RT_iterator>
make_cartesian_const_iterator_begin(RT_iterator begin, RT_iterator w)
{
  return Cartesian_const_iterator_d<RT_iterator>(begin, w);
}

template < typename RT_iterator > inline
Cartesian_const_iterator_d<RT_iterator>
make_cartesian_const_iterator_end(RT_iterator w)
{
  return Cartesian_const_iterator_d<RT_iterator>(w, w);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_CONST_ITERATOR_D_H
