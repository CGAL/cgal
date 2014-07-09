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

  const I& begin() const
  {
    return this->first;
  }

  const I& end() const
  {
    return this->second;
  }
};

} // namespace CGAL

#endif // CGAL_ITERATOR_RANGE_H
