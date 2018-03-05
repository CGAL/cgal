// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Maxime Gimeno

#ifndef GEOMETRY_CONTAINER_H
#define GEOMETRY_CONTAINER_H
namespace CGAL{

//!
//! \brief The Geometry_container struct is a wrapper 
//! that provides an easy way to use the read and write functions
//! of `Polygon_IO`.
//! 
//! \tparam Range is a model of `RandomAccessRange`
//! \tparam TAG is the tag corresponding to the wkt structure you want
//! to be read or written with this `CGAL::Geometry_container`
//!
template <typename Range, typename TAG>
struct Geometry_container{
  typedef std::pair<Range, TAG> type;
  typedef typename Range::difference_type difference_type; 
  typedef typename Range::iterator iterator;
  typedef typename Range::const_iterator const_iterator;
  typedef typename Range::reverse_iterator reverse_iterator;
  typedef typename Range::const_reverse_iterator const_reverse_iterator;
  typedef typename Range::size_type size_type;
  typedef typename Range::value_type value_type;

  Range range;
  //!
  //! \brief default constructor. 
  //!
  //! It will create an internal  default `Range`.
  //!
  Geometry_container() {}
  //!
  //! \brief copy constructor.
  //! 
  //!It will copy `range` into the internal `Range`.
  //!
  Geometry_container(Range& range)
    :range(range){}
  
  iterator begin()
  { return range.begin(); }
  
  iterator end()
  { return range.end(); }
  
  reverse_iterator rbegin()
  { return range.rbegin(); }
  
  reverse_iterator rend()
  { return range.rend(); }
  
  const_iterator begin()const
  { return range.begin(); }
  
  const_iterator end()const
  { return range.end(); }
  
  const_reverse_iterator rbegin()const
  { return range.rbegin(); }
  
  const_reverse_iterator rend()const
  { return range.rend(); }
  
  void clear(){ range.clear(); }
  
  template<typename size_t>
  void resize(size_t n){ range.resize(n); }
  
  template<typename T>
  void push_back(const T& t){ range.push_back(t); }
  
  size_type size() { return range.size(); }
  
  value_type operator[](size_type i)
  { return range[i]; }
};//end Geometry_container
}//end CGAL

namespace boost{

template< class T, typename TAG >
struct range_iterator<CGAL::Geometry_container<T, TAG> >
{ typedef typename T::iterator type; };

template< class T, typename TAG >
struct range_iterator<const CGAL::Geometry_container<T, TAG> >
{ typedef typename T::const_iterator type; };

template< class T, typename TAG >
struct range_mutable_iterator<CGAL::Geometry_container<T, TAG> >
{ typedef typename range_mutable_iterator<T>::type type; };

}//end boost
#endif // GEOMETRY_CONTAINER_H
