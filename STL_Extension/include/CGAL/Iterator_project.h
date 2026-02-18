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

#include <CGAL/config.h>
#include <CGAL/type_traits.h>

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

template <class I, class Fct>
using Iterator_project_reference = decltype(std::declval<Fct>()(*std::declval<I>()));

template <class I, class Fct>
using Iterator_project_value_type = CGAL::cpp20::remove_cvref_t<Iterator_project_reference<I, Fct>>;

// keep 4 dummy template parameters around for backwards compatibility
template < class I, class Fct,
           class D1 = int, class D2 = int, class D3 = int, class D4 = int >
class Iterator_project
  : public boost::stl_interfaces::v1::proxy_iterator_interface<
               Iterator_project<I, Fct, D1, D2, D3, D4>,
               typename std::iterator_traits<I>::iterator_category,
               Iterator_project_value_type<I, Fct>,
               Iterator_project_reference<I, Fct>
               >
{
public:
  using reference = Iterator_project_reference<I, Fct>;
  using value_type = Iterator_project_value_type<I, Fct>;
protected:
  using Base = boost::stl_interfaces::v1::proxy_iterator_interface<
      Iterator_project<I, Fct, D1, D2, D3, D4>,
      typename std::iterator_traits<I>::iterator_category,
      value_type,
      reference>;

  I        nt;    // The internal iterator.

  friend boost::stl_interfaces::access;

  I& base_reference() { return nt; }
  const I& base_reference() const { return nt; }

public:
  using Self = Iterator_project<I,Fct,D1,D2,D3,D4>;
  using Iterator = I; // base iterator
  using pointer = typename Base::pointer;
  using difference_type = typename Base::difference_type;

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

  reference operator*() const {
    Fct fct;
    return fct(*nt);
  }

  // AF: Added in this PR
  pointer ptr() const {
    return this->operator->();
  }
  // OPERATIONS Random Access Category
  // ---------------------------------

  template <typename T = I>
  auto operator+=( difference_type n) -> decltype(std::declval<T&>() += n, std::declval<Self&>()) {
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

#if defined(BOOST_MSVC)
  difference_type operator- (const Self& i) const {
    return nt - i.nt;
  }
#else

template <typename It>
  std::enable_if_t<std::is_convertible_v<typename It::iterator_category, std::random_access_iterator_tag>
                   && std::is_convertible_v<const It&, const Self&>,
                   difference_type>
  operator- (const It& i) const {
    return nt - i.nt;
  }
#endif

  reference operator[]( difference_type n) const {
    Self tmp = *this;
    tmp += n;
    return tmp.operator*();
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
