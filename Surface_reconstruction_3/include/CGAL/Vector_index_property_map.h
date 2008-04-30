// Copyright (c) 2007-2008  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Laurent Saboret

#ifndef CGAL_VECTOR_INDEX_PROPERTY_MAP_H
#define CGAL_VECTOR_INDEX_PROPERTY_MAP_H

#include <boost/property_map.hpp>

#include <vector>

CGAL_BEGIN_NAMESPACE


/// Type of the "vertex_index" property map
/// of a std::vector object.
template <class T>
class Vector_index_property_map 
{
public:
    // Property maps required types
    typedef boost::readable_property_map_tag          category;
    typedef unsigned int                              value_type;
    typedef value_type                                reference;
    typedef typename std::vector<T>::const_iterator   key_type;

    Vector_index_property_map(const std::vector<T>& data) 
    : m_data(data) 
    {}

    /// Free function to access the map elements.
    friend inline
    reference get(const Vector_index_property_map& map, key_type p)
    {
      return std::distance(map.m_data.begin(), p);
    }
    
private:
  const std::vector<T>& m_data;
};

/// Free function to get the "vertex_index" property map
/// of a std::vector object.
template <class T>
inline
Vector_index_property_map<T> 
get(boost::vertex_index_t, const std::vector<T>& data) 
{
  Vector_index_property_map<T> aMap(data);
  return aMap;
}


CGAL_END_NAMESPACE

#endif //CGAL_VECTOR_INDEX_PROPERTY_MAP_H

