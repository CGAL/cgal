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

#ifndef CGAL_RANDOM_ACCESS_ITERATOR_INDEX_PMAP_H
#define CGAL_RANDOM_ACCESS_ITERATOR_INDEX_PMAP_H

#include <boost/property_map.hpp>
#include <boost/graph/properties.hpp>

#include <vector>
#include <deque>

CGAL_BEGIN_NAMESPACE


/// Type of the "vertex_index" property map
/// of a random access container (typically vector and deque).
template <class RandomAccessContainer>
class Random_access_container_index_pmap 
{
public:
    // Property maps required types
    typedef boost::readable_property_map_tag                category;
    typedef unsigned int                                    value_type;
    typedef value_type                                      reference;
    typedef typename RandomAccessContainer::const_iterator  key_type;

    Random_access_container_index_pmap(const RandomAccessContainer& data) 
    : m_data(data) 
    {}

    /// Free function to access the map elements.
    friend inline
    reference get(const Random_access_container_index_pmap& map, key_type p)
    {
      // Safety: the next line is dentical to std::distance(map.m_data.begin(), p)
      // but will fail to compile for non random access containers.
      return p - map.m_data.begin();
    }
    
private:
  const RandomAccessContainer& m_data;
};


/// Free function to get the "vertex_index" property map
/// of a std::vector object.
template <class T>
inline
Random_access_container_index_pmap< std::vector<T> > 
get(boost::vertex_index_t, const std::vector<T>& data) 
{
  Random_access_container_index_pmap< std::vector<T> > aMap(data);
  return aMap;
}

/// Free function to get the "vertex_index" property map
/// of a std::deque object.
template <class T>
inline
Random_access_container_index_pmap< std::deque<T> > 
get(boost::vertex_index_t, const std::deque<T>& data) 
{
  Random_access_container_index_pmap< std::deque<T> > aMap(data);
  return aMap;
}


CGAL_END_NAMESPACE

#endif //CGAL_RANDOM_ACCESS_ITERATOR_INDEX_PMAP_H

