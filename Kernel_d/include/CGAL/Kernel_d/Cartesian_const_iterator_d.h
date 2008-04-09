// Copyright (c) 2000,2001,2008  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Kernel_d/include/CGAL/Kernel_d/Tuple_d.h $
// $Id: Tuple_d.h 42814 2008-04-09 16:07:00Z spion $
//
// Author(s)     : Michael Seel, Sylvain Pion

#ifndef CGAL_CARTESIAN_CONST_ITERATOR_D_H
#define CGAL_CARTESIAN_CONST_ITERATOR_D_H

#include <CGAL/basic.h>
#include <CGAL/Quotient.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE

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

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONST_ITERATOR_D_H
