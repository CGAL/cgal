// Copyright (c) 2007 
// GeometryFactory (France),
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel). All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion
//                 Fernando Cacciola <fernando.cacciola@geometryfactory.com> 

#ifndef CGAL_ITERATOR_TRANSFORM_H
#define CGAL_ITERATOR_TRANSFORM_H 1

#include <CGAL/Iterator_project.h>

namespace CGAL {

template < class I, class Fct>
class Iterator_transform {
protected:
  I        nt;    // The internal iterator.
public:
  typedef Iterator_transform<I,Fct> Self;
  typedef I                                   Iterator; // base iterator
  typedef std::iterator_traits<I>             traits;
  typedef typename traits::difference_type    difference_type;
  typedef typename traits::iterator_category  iterator_category;
  typedef typename traits::value_type         base_value_type;
  typedef typename traits::pointer            base_pointer;
  typedef typename traits::reference          base_reference;

  typedef typename Fct::argument_type         argument_type;
  typedef typename Fct::result_type           value_type;


  // This iterator returns rvalues by design (allowing the conversion function to return new objects)
  typedef value_type reference;

  // Use I_TYPE_MATCH_IF to find correct pointer type.
  typedef I_TYPE_MATCH_IF< base_pointer, const base_value_type *,
    const value_type *, value_type *> Match2;
  typedef typename Match2::Result             pointer;

  // CREATION
  // --------

  Iterator_transform() {}
  Iterator_transform( I j) : nt(j) {}

  // make two iterators assignable if the underlying iterators are
  template <class I2>
  Iterator_transform( const Iterator_transform<I2,Fct>& i2)
  : nt( i2.current_iterator()) {}

  template <class I2>
  Self& operator= ( const Iterator_transform<I2,Fct>& i2) {
    nt = i2.current_iterator();
    return *this;
  }

  // OPERATIONS Forward Category
  // ---------------------------

  Iterator  current_iterator() const { return nt;}
  bool      operator==( const Self& i) const { return ( nt == i.nt); }
  bool      operator!=( const Self& i) const { return !(*this == i); }
  
  struct Proxy
  {
      Proxy(const reference r) : ref(r) {}
      reference ref;
      pointer operator->() { return &ref; }
  };

  Proxy operator->() const
  {
      return Proxy(Fct()(*nt));
  }
  
  reference operator* () const 
  {
    Fct fct;
    return fct(*nt);
  }
  
  Self&     operator++() {
    ++nt;
    return *this;
  }
  Self      operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  // OPERATIONS Bidirectional Category
  // ---------------------------------

  Self& operator--() {
    --nt;
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }

  // OPERATIONS Random Access Category
  // ---------------------------------

  Self& operator+=( difference_type n) {
    nt += n;
    return *this;
  }
  Self  operator+( difference_type n) const {
    Self tmp = *this;
    return tmp += n;
  }
  Self& operator-=( difference_type n) {
    return operator+=( -n);
  }
  Self  operator-( difference_type n) const {
    Self tmp = *this;
    return tmp += -n;
  }
  difference_type  operator-( const Self& i) const { return nt - i.nt; }
  reference  operator[]( difference_type n) const {
    Self tmp = *this;
    tmp += n;
    return tmp.operator*();
  }
  bool operator< ( const Self& i) const { return ( nt < i.nt); }
  bool operator> ( const Self& i) const { return i < *this; }
  bool operator<=( const Self& i) const { return !(i < *this); }
  bool operator>=( const Self& i) const { return !(*this < i); }
};

template < class Dist, class Fct, class I>
inline
Iterator_transform<I,Fct>
operator+( Dist n, Iterator_transform<I,Fct> i) {
  return i += n;
}

} //namespace CGAL
#endif // CGAL_Iterator_transform_H //
// EOF //
