// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Ophir Setter <ophirset@post.tau.ac.il>

/*! \file
 * This file handles the color property_map for running certain boost
 * algorithms (e.g., BFS and DFS) on data structures that use handles to
 * store incident relations. The map values, that are the colors, are stored
 * at the handled object, and directly accessed through the handle.
 */

#ifndef CGAL_OBJECT_COLOR_MAP_H
#define CGAL_OBJECT_COLOR_MAP_H

#include <boost/property_map.hpp>
#include <boost/graph/properties.hpp>

CGAL_BEGIN_NAMESPACE

/*! \class
 * This class can be used as the color property-map when using certain boost
 * algorithms (e.g., BFS and DFS) on data structures that use handles to
 * store incident relations. Given a handle, this property map stores and
 * retrieves its color property in constant time. It assumes that the property
 * can be stored and retrieved through the handle with the data() and
 * set_data() methods respectively.
 * The class models the ReadWritePropertyMapConcept concept, where the
 * map property key is a handle, and its value is a default_color_type.
 *
 * A private case of this is a face that is created using the Arr_extended_dcel
 * object.
 */
template <class Handle>
class Object_color_map
{
public:
  typedef boost::read_write_property_map_tag    category;
  typedef boost::default_color_type             value_type;
  typedef boost::default_color_type             reference;
  typedef Handle                                key_type;

  // get and set function to get and set the color.
  value_type get(key_type handle) const { return (*handle).data(); }
  void put(key_type & handle, const value_type & color) 
  { 
    (*handle).set_data(color); 
  }
};

/*!
 * Specialize the get property_map function, so that BOOST is able to use the
 * above class.
 * @param col_map The color map to take the color from.
 * @param k The handle to the colored object.
 * @return The color of the handled object.
 */
template<class Handle>
boost::default_color_type get(Object_color_map<Handle> col_map, Handle k) 
{ 
  return col_map.get(k);
}

/*!
 * Specialize the set property_map function, so that BOOST is able to use the
 * above class.
 * @param col_map The color map to take the color from.
 * @param val The color the coloer object is to be set with.
 */
template<class Handle>
void put(Object_color_map<Handle> col_map, Handle k,
         boost::default_color_type val)
{ 
  col_map.put( k, val );
}

CGAL_END_NAMESPACE

#endif
