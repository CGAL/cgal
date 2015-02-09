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

#endif // CGAL_ITERATOR_RANGE_H
