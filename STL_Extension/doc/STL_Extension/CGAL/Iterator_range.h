// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_ITERATOR_RANGE_H
#define CGAL_ITERATOR_RANGE_H

#include <CGAL/tuple.h>
#include <utility>

namespace CGAL {

  /*!
    \ingroup PkgSTLExtensionRef
    `CGAL::Iterator_range` encapsulates two iterators so they fulfill the `ForwardRange` concept.
    The class is essentially a clone of <A href="https://www.boost.org/doc/libs/1_55_0/libs/range/doc/html/range/reference/utilities/iterator_range.html">`boost::iterator_range`</A>,
    and it additionally is derived from `std::pair`, so that one can apply `boost::tie`.
  */
  template <typename I>
  class Iterator_range
    : public std::pair<I,I>{

    typedef std::pair<I,I> Base;

  public:

    typedef I iterator;
    typedef I const_iterator;

    Iterator_range(I b, I e)
      : Base(b,e)
    {}

    Iterator_range(const std::pair<I,I>& ip)
      : Base(ip)
    {}


  I begin() const
  {
    return this->first;
  }

  I end() const
  {
    return this->second;
  }

  /// returns `std::distance(begin(), end())`
  std::size_t
  size() const
  {
    return static_cast<std::size_t>(std::distance(begin(), end()));
  }
};

  template <typename T>
  Iterator_range<T>
  make_range(const T& b, const T&e)
  {
    return Iterator_range<T>(b,e);
  }

  template<typename T>
  inline T range_begin( Iterator_range<T> & x )
  {
    return x.first;
  }

  template<typename T>
  inline T range_end( Iterator_range<T> & x )
  {
    return x.second;
  }

  template<typename T>
  inline T range_begin(const Iterator_range<T>& x )
  {
    return x.first;
  }

  template<typename T>
  inline T range_end(const Iterator_range<T>& x )
  {
    return x.second;
  }
} // namespace CGAL

#ifndef DOXYGEN_RUNNING

namespace boost {

  template <typename X>
  struct range_iterator;

  template <typename X>
  struct range_mutable_iterator;

  template <typename X>
  struct range_const_iterator;

  template <typename T>
  struct range_iterator<CGAL::Iterator_range<T> >
  {
    typedef T type;
  };


  template<typename T>
  struct range_mutable_iterator< CGAL::Iterator_range<T> >
  {
    typedef T type;
  };

  template<typename T>
  struct range_const_iterator< CGAL::Iterator_range<T> >
  {
    typedef T type;
  };
} // namespace boost

#endif // DOXYGEN_RUNNING

#endif // CGAL_ITERATOR_RANGE_H
