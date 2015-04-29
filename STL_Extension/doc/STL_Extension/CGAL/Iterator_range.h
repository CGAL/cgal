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
// 
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_ITERATOR_RANGE_H
#define CGAL_ITERATOR_RANGE_H

#include <CGAL/tuple.h>
#include <utility>

namespace CGAL {

  /*!
    \ingroup PkgStlExtension
    `CGAL::Iterator_range` encapsulates two iterators so they fulfill the `ForwardRange` concept. 
    The class is essentially a clone of <A href="http://www.boost.org/doc/libs/1_55_0/libs/range/doc/html/range/reference/utilities/iterator_range.html">`boost::iterator_range`</A>,
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
  typename std::iterator_traits<I>::difference_type
  size() const
  {
    return std::distance(begin(), end());
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
}
#endif // CGAL_ITERATOR_RANGE_H
