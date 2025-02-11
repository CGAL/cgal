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
#include <boost/foreach.hpp>

namespace CGAL {

/*!
  \ingroup PkgSTLExtensionRef
  `CGAL::Iterator_range` encapsulates two iterators so they fulfill the `ForwardRange` concept.
  The class is essentially a clone of <A
  href="https://www.boost.org/doc/libs/1_55_0/libs/range/doc/html/range/reference/utilities/iterator_range.html">`boost::iterator_range`</A>,
  and it additionally is derived from `std::pair`, so that one can apply `boost::tie`.
*/
template <typename I>
class Iterator_range : public std::pair<I, I>
{
  typedef std::pair<I, I> Base;

public:
  typedef I iterator;
  typedef I const_iterator;

  Iterator_range() = default;
  Iterator_range(I b, I e)
      : Base(b, e) {}
  Iterator_range(const std::pair<I, I>& ip)
      : Base(ip) {}

  I begin() const { return this->first; }
  I end() const { return this->second; }

  /// returns `std::distance(begin(), end())`
  std::size_t size() const { return static_cast<std::size_t>(std::distance(begin(), end())); }

  /// returns `std::distance(begin(), end())==0`
  bool empty() const { return begin() == end(); }

  operator std::tuple<I&, I&>()
  {
    return std::tuple<I&, I&>{this->first, this->second};
  }

  operator std::tuple<const I&, const I&>() const
  {
    return std::tuple<const I&, const I&>{this->first, this->second};
  }

  template <template<class...> class Container>
  auto to() const
  {
    using V = std::remove_cv_t<std::remove_reference_t<decltype(*begin())>>;
    return Container<V>(begin(), end());
  }
};

template <typename T>
Iterator_range<T>
make_range(const T& b, const T&e)
{
  return Iterator_range<T>(b,e);
}

template <typename T>
Iterator_range<T>
make_range(const std::pair<T,T>& p)
{
  return Iterator_range<T>(p.first,p.second);
}


} // namespace CGAL

namespace boost::foreach {
  template<typename T>
  struct is_lightweight_proxy<CGAL::Iterator_range<T>> : boost::mpl::true_ {};
}

#if CGAL_CXX20
#  include <ranges>

  template<typename I>
  inline constexpr bool std::ranges::enable_borrowed_range<CGAL::Iterator_range<I>> = true;

#endif // C++20

#endif // CGAL_ITERATOR_RANGE_H
