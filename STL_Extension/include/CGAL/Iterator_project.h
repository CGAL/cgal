// Copyright (c) 2003
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion

#ifndef CGAL_ITERATOR_PROJECT_H
#define CGAL_ITERATOR_PROJECT_H 1

#include <iterator>

namespace CGAL {

// Relies on iterator traits. Quite simplified compared to earlier version.

// The pointer type and the reference type in the Iterator_project
// are based on the value type from the projector, but the base iterator
// determines whether they are const or mutable. The following template
// class and its partial specialization helps creating the derived types.

// If T === T1 return R1 else return R2
template <class T, class T1, class R1, class R2>
struct I_TYPE_MATCH_IF { typedef R2 Result; };  // else clause

template <class T, class R1, class R2>
struct I_TYPE_MATCH_IF<T,T,R1,R2> { typedef R1 Result; }; // then clause

// keep 4 dummy template parameters around for backwards compatibility
template < class I, class Fct,
           class D1 = int, class D2 = int, class D3 = int, class D4 = int >
class Iterator_project {
protected:
  I        nt;    // The internal iterator.
public:
  typedef Iterator_project<I,Fct,D1,D2,D3,D4> Self;
  typedef I                                   Iterator; // base iterator
  typedef std::iterator_traits<I>             traits;
  typedef typename traits::difference_type    difference_type;
  typedef typename traits::iterator_category  iterator_category;
  typedef typename traits::value_type         base_value_type;
  typedef typename traits::pointer            base_pointer;
  typedef typename traits::reference          base_reference;

  typedef typename Fct::argument_type         argument_type;
  typedef typename Fct::result_type           value_type;

  // Use I_TYPE_MATCH_IF to find correct pointer and reference type.

  typedef I_TYPE_MATCH_IF< base_reference, const base_value_type &,
    const value_type &, value_type &> Match1;
  typedef typename Match1::Result             reference;

  typedef I_TYPE_MATCH_IF< base_pointer, const base_value_type *,
    const value_type *, value_type *> Match2;
  typedef typename Match2::Result             pointer;

  // CREATION
  // --------

  Iterator_project() {}
  Iterator_project( I j) : nt(j) {}

  // make two iterators assignable if the underlying iterators are
  template <class I2, class Q1, class Q2, class Q3, class Q4>
  Iterator_project( const Iterator_project<I2,Fct,Q1,Q2,Q3,Q4>& i2)
  : nt( i2.current_iterator()) {}

  template <class I2, class Q1, class Q2, class Q3, class Q4>
  Self& operator= ( const Iterator_project<I2,Fct,Q1,Q2,Q3,Q4>& i2) {
    nt = i2.current_iterator();
    return *this;
  }

  // OPERATIONS Forward Category
  // ---------------------------

  Iterator  current_iterator() const { return nt;}
  pointer   ptr() const {
    Fct fct;
    return &(fct(*nt));
  }
  bool      operator==( const Self& i) const { return ( nt == i.nt); }
  bool      operator!=( const Self& i) const { return !(*this == i); }
  reference operator* ()               const { return *ptr(); }
  pointer   operator->()               const { return ptr(); }
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

template < class Dist, class Fct, class I,
           class D1, class D2, class D3, class D4>
inline
Iterator_project<I,Fct,D1,D2,D3,D4>
operator+( Dist n, Iterator_project<I,Fct,D1,D2,D3,D4> i) {
  return i += n;
}

} //namespace CGAL
#endif // CGAL_ITERATOR_PROJECT_H //
// EOF //
