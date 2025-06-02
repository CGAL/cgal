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
#include <boost/stl_interfaces/iterator_interface.hpp>

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
class Iterator_project
  : public boost::stl_interfaces::proxy_iterator_interface<Iterator_project<I, Fct, D1, D2, D3, D4>,
                                                           typename std::iterator_traits<I>::iterator_category,
                                                           typename Fct::result_type>
{
  using Base = boost::stl_interfaces::proxy_iterator_interface<
      Iterator_project<I, Fct, D1, D2, D3, D4>,
      typename std::iterator_traits<I>::iterator_category,
      typename Fct::result_type>;

protected:
  I        nt;    // The internal iterator.

  friend boost::stl_interfaces::access;

  I& base_reference() { return nt; }
  const I& base_reference() const { return nt; }
public:
  using Self = Iterator_project<I,Fct,D1,D2,D3,D4>;
  using Iterator = I; // base iterator
  using value_type = typename Fct::result_type;
  using pointer = typename Base::pointer;

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

  value_type operator*() const {
    Fct fct;
    return fct(*nt);
  }

  pointer ptr() const {
    return this->operator->();
  }
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
