// Copyright (c) 2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_COORDINATE_ITERATOR_2_H
#define CGAL_CARTESIAN_COORDINATE_ITERATOR_2_H

#include <cstddef>
#include <iterator>

CGAL_BEGIN_NAMESPACE


template <class K>
class Cartesian_coordinate_iterator_2 
{

protected:
  typedef typename K::Point_2 P;
  const P* point;
  int index;
  typedef Cartesian_coordinate_iterator_2<K> Self;

public: 

  typedef typename K::FT FT;
  typedef P Point;

  typedef std::random_access_iterator_tag iterator_category;
  typedef FT                              value_type;
  typedef std::ptrdiff_t                  difference_type;
  typedef const value_type&               reference;
  typedef const value_type*               pointer;

  Cartesian_coordinate_iterator_2(const Point *const p = NULL, 
				  int _index = 0) 
    : point(p), index(_index) 
  {}


  const FT 
  operator*() const {
    return point->cartesian(index); 
  }
  
  Self&  operator++() {
    index++; 
    return *this;
  }

  Self&  
  operator--() {
    index--; 
    return *this;
  }

  Self 
  operator++(int) { 
    Self tmp(*this);
    ++(*this); 
    return tmp; 
  }

  Self
  operator--(int) {
    Self tmp(*this);
    --(*this);
    return tmp;
  }
  
  Self& 
  operator+=(difference_type i) { 
    index+=i;
    return *this; 
  }

  Self& 
  operator-=(difference_type i) { 
    index -= i; 
    return *this; 
  }

  Self 
  operator+(difference_type i) const {
    Self tmp=*this; 
    return tmp += i; 
  }

  Self operator-(difference_type i) const {
    Self tmp=*this; 
    return tmp -= i; 
  }

  difference_type 
  operator-(const Self& x) const {
    CGAL_kernel_assertion(point == x.point);
    return index - x.index; 
  }

  reference operator[](difference_type i) const { 
    return *(*this + i); 
  }

  bool operator==(const Self& x) const {
    return (point == x.point)&& (index == x.index) ; 
  }
  
  bool operator!=(const Self& x) const { 
    return ! (*this==x); 
  }

  bool operator<(const Self& x) const 
  { 
    return (x - *this) > 0; }

};



CGAL_END_NAMESPACE
#endif // CGAL_CARTESIAN_COORDINATE_ITERATOR_2_H
