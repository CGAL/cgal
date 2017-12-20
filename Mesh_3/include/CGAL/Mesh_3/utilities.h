// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_3_UTILITIES_H
#define CGAL_MESH_3_UTILITIES_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Has_timestamp.h>
#include <iterator>
#include <string>
#include <sstream>

namespace CGAL {

namespace Mesh_3 {
namespace internal {
  
struct Debug_messages_tools {
  template <typename Vertex_handle>
  static std::string disp_vert(Vertex_handle v, Tag_true) {
    std::stringstream ss;
    ss.precision(17);
    ss << (void*)(&*v) << "[ts=" << v->time_stamp() << "]"
       << "(" << v->point() <<")";
    return ss.str();
  }

  template <typename Vertex_handle>
  static std::string disp_vert(Vertex_handle v, Tag_false) {
    std::stringstream ss;
    ss.precision(17);
    ss << (void*)(&*v) << "(" << v->point() <<")";
    return ss.str();
  }

  template <typename Vertex_handle>
  static std::string disp_vert(Vertex_handle v)
  {
    typedef typename std::iterator_traits<Vertex_handle>::value_type Vertex;
    return disp_vert(v, CGAL::internal::Has_timestamp<Vertex>());
  }
};

/**
 * @class First_of
 * Function object which returns the first element of a pair
 */
template <typename Pair>
struct First_of :
  public CGAL::unary_function<Pair, const typename Pair::first_type&>
{
  typedef CGAL::unary_function<Pair, const typename Pair::first_type&> Base;
  typedef typename Base::result_type                                  result_type;
  typedef typename Base::argument_type                                argument_type;
  
  result_type operator()(const argument_type& p) const { return p.first; }
}; // end class First_of
  

/**
 * @class Ordered_pair
 * Stores two elements in an ordered manner, i.e. first() < second()
 */
template <typename T>
class Ordered_pair
{
public:
  Ordered_pair(const T& t1, const T& t2)
  : data_(t1,t2)
  {
    if ( ! (t1 < t2) ) 
    {
      data_.second = t1;
      data_.first = t2;
    }
  }
  
  const T& first() const { return data_.first; }
  const T& second() const { return data_.second; }
  
  bool operator<(const Ordered_pair& rhs) const { return data_ < rhs.data_; }
  
private:
  std::pair<T,T> data_;
};


/**
 * @class Iterator_not_in_complex
 * @brief A class to filter elements which do not belong to the complex
 */
template < typename C3T3 >
class Iterator_not_in_complex
{
  const C3T3& c3t3_;
public:
  Iterator_not_in_complex(const C3T3& c3t3) : c3t3_(c3t3) { }
  
  template <typename Iterator>
  bool operator()(Iterator it) const { return ! c3t3_.is_in_complex(*it); }
}; // end class Iterator_not_in_complex


} // end namespace internal  
} // end namespace Mesh_3


} //namespace CGAL

#endif // CGAL_MESH_3_UTILITIES_H
