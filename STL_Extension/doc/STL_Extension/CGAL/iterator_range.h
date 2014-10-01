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
    `CGAL::iterator_range` is a...
  */
  template <typename I>
  class iterator_range
    : public std::pair<I,I>{
    
    typedef std::pair<I,I> Base;

  public:

    typedef I iterator;
    typedef const I const_iterator;

    iterator_range(I b, I e)
      : Base(b,e)
    {}

    iterator_range(const std::pair<I,I>& ip)
      : Base(ip)
    {}

    operator Base() const
    {
      return std::make_pair(begin(),end());
    }

  const I& begin() const
  {
    return this->first;
  }

  const I& end() const
  {
    return this->second;
  }
};

  template <typename T>
  iterator_range<T>
  make_range(const T& b, const T&e)
  {
    return iterator_range<T>(b,e);
  }

  template<typename T>
  inline T range_begin( iterator_range<T> & x )
  {
    return x.first;
  }
  
  template<typename T>
  inline T range_end( iterator_range<T> & x )
  {
    return x.second;
  }
  
  template<typename T>
  inline T range_begin(const iterator_range<T>& x )
  {
    return x.first;
  }
  
  template<typename T>
  inline T range_end(const iterator_range<T>& x )
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
  struct range_iterator<CGAL::iterator_range<T> >
  {
    typedef T type;
  };


  template<typename T>
  struct range_mutable_iterator< CGAL::iterator_range<T> >
  {
    typedef T type;
  };
  
  template<typename T>
  struct range_const_iterator< CGAL::iterator_range<T> >
  {
    typedef T type;
  };
}
#endif // CGAL_ITERATOR_RANGE_H
