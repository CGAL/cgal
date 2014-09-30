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
// Author(s)     : Andreas Fabri
 
#ifndef CGAL_RANGE_H
#define CGAL_RANGE_H

#include <utility>

namespace CGAL {

// Range is a ... 

  template <typename T>
  struct Range : public std::pair<T,T>
  {
    typedef std::pair<T,T> Base;

    typedef T iterator;

    Range(const T&b, const T& e)
      : Base(b,e)
    {}

    const T& begin() const 
    {
      return first;
    }  

    const T& end() const 
    {
      return second;
    }
  };

  template <typename T>
  Range<T>
  make_range(const T& b, const T&e)
  {
    return Range<T>(b,e);
  }

  template<typename T>
  inline T range_begin( Range<T> & x )
  {
    return x.first;
  }
  
  template<typename T>
  inline T range_end( Range<T> & x )
  {
    return x.second;
  }
  
  template<typename T>
  inline T range_begin(const Range<T>& x )
  {
    return x.first;
  }
  
  template<typename T>
  inline T range_end(const Range<T>& x )
  {
    return x.second;
  }  
}

namespace boost {
  
  template <typename X> 
  struct range_iterator;

  template <typename X> 
  struct range_mutable_iterator;

  template <typename X> 
  struct range_const_iterator;

  template <typename T> 
  struct range_iterator<CGAL::Range<T> >
  {
    typedef T type;
  };


  template<typename T>
  struct range_mutable_iterator< CGAL::Range<T> >
  {
    typedef T type;
  };
  
  template<typename T>
  struct range_const_iterator< CGAL::Range<T> >
  {
    typedef T type;
  };
  
  /*  
  template <typename T>
  T
  begin(const CGAL::Range<T>& r)
  {
    return r.first;
  }
  
  template <typename T>
  T
  end(const CGAL::Range<T>& r)
  {
    return r.second;
  }
  */
}

#endif // CGAL_RANGE_H
