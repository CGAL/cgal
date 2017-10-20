// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_ITERATOR_RANGE_H
#define CGAL_ITERATOR_RANGE_H

#include <CGAL/tuple.h>
#include <utility>
#include <boost/foreach.hpp>

namespace CGAL {

  /*!
\ingroup PkgStlExtension
  /// `CGAL::Iterator_range` is a...
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


    // Iterator_range(const Iterator_range& ip)
    //   : Base(ip)
    // {}

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

  /// returns `std::distance(begin(), end())==0`
  bool empty() const
  {
    return begin()==end();
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

// At global scope...

  template<typename T>
inline boost::mpl::true_ *
  boost_foreach_is_lightweight_proxy( CGAL::Iterator_range<T> *&, boost::foreach::tag )
{
    return 0;
}
namespace boost { namespace foreach
{
    template<typename T>
    struct is_lightweight_proxy< CGAL::Iterator_range<T> >
      : mpl::true_
    {
    };
}}
#endif // CGAL_ITERATOR_RANGE_H
