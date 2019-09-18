// Copyright (c) 2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_VERTEX_INDEX_MAP_H
#define CGAL_ARR_VERTEX_INDEX_MAP_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>


/*! \file
 * Definition of the Arr_vertex_index_map<Arrangement> class.
 */
#include <CGAL/Arr_observer.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/property_map.h>

#include <boost/graph/properties.hpp>

namespace CGAL {

/*! \class
 * An auxiliary class that automatically maintains a mapping of the
 * arrangement vertices to the indices 0, ..., (n -1), where n is the number
 * of vertices in the arrangement.
 */
template <class Arrangement_>
class Arr_vertex_index_map : public Arr_observer<Arrangement_>
{
public:

  typedef Arrangement_                            Arrangement_2;
  typedef typename Arrangement_2::Vertex_handle   Vertex_handle;

  // Boost property type definitions:
  typedef boost::readable_property_map_tag        category;
  typedef unsigned int                            value_type;
  typedef value_type                              reference;
  typedef Vertex_handle                           key_type;

private:

  typedef Arr_vertex_index_map<Arrangement_2>     Self;
  typedef Arr_observer<Arrangement_2>             Base;

  typedef Unique_hash_map<Vertex_handle, unsigned int>     Index_map;

  // Data members:
  unsigned int                n_vertices;  // The current number of vertices.
  Index_map                   index_map;   // Mapping vertices to indices.
  std::vector<Vertex_handle>  rev_map;     // Mapping indices to vertices.

  enum {MIN_REV_MAP_SIZE = 32};

public:

  /*! Default constructor. */
  Arr_vertex_index_map () :
    Base (),
    n_vertices (0),
    rev_map (MIN_REV_MAP_SIZE)
  {}

  /*! Constructor with an associated arrangement. */
  Arr_vertex_index_map (const Arrangement_2& arr) :
    Base (const_cast<Arrangement_2&> (arr))
  {
    _init();
  }

  /*! Copy constructor. */
  Arr_vertex_index_map (const Self& other) :
    Base (const_cast<Arrangement_2&> (*(other.arrangement())))
  {
    _init();
  }

  /*! Assignment operator. */
  Self& operator= (const Self& other)
  {
    if (this == &other)
      return (*this);

    this->detach();
    this->attach (const_cast<Arrangement_2&> (*(other.arrangement())));

    return (*this);
  }

  /*!
   * Get the index of a given vertex.
   * \param v A handle to the vertex.
   * \pre v is a valid vertex in the arrangement.
   */
  unsigned int operator[] (Vertex_handle v) const
  {
    return index_map[v];
  }

  /*!
   * Get the vertex given its index.
   * \param i The index of the vertex.
   * \pre i is less than the number of vertices in the graph.
   */
  Vertex_handle vertex (const int i) const
  {
    CGAL_precondition (i < n_vertices);

    return rev_map[i];
  }

  /// \name Notification functions, to keep the mapping up-to-date.
  //@{

  /*!
   * Update the mapping after the arrangement has been assigned with another
   * arrangement.
   */
  virtual void after_assign ()
  {
    _init();
  }

  /*!
   * Update the mapping after the arrangement is cleared.
   */
  virtual void after_clear ()
  {
    _init();
  }

  /*!
   * Update the mapping after attaching to a new arrangement.
   */
  virtual void after_attach ()
  {
    _init();
  }

  /*!
   * Update the mapping after detaching the arrangement.
   */
  virtual void after_detach ()
  {
    n_vertices = 0;
    index_map.clear();
  }

  /*!
   * Update the mapping after the creation of a new vertex.
   * \param v A handle to the created vertex.
   */
  virtual void after_create_vertex (Vertex_handle v)
  {
    // Update the number of vertices.
    n_vertices++;

    // If necessary, allocate memory for the reverse mapping.
    if (rev_map.size() < n_vertices)
      rev_map.resize (2 * n_vertices);

    // Update the mapping of the newly created vertex.
    index_map[v] = n_vertices - 1;
    rev_map[n_vertices - 1] = v;
  }

  /*!
   * Update the mapping after the creation of a new boundary vertex.
   * \param v A handle to the created vertex.
   */
  virtual void after_create_boundary_vertex (Vertex_handle v)
  {
    // Update the number of vertices.
    n_vertices++;

    // If necessary, allocate memory for the reverse mapping.
    if (rev_map.size() < n_vertices)
      rev_map.resize (2 * n_vertices);

    // Update the mapping of the newly created vertex.
    index_map[v] = n_vertices - 1;
    rev_map[n_vertices - 1] = v;
  }

  /*!
   * Update the mapping before the removal of a vertex.
   * \param v A handle to the vertex to be removed.
   */
  virtual void before_remove_vertex (Vertex_handle v)
  {
    // Update the number of vertices.
    n_vertices--;
    
    // Reduce memory consumption in case the number of vertices has
    // drastically decreased.
    if (2*n_vertices+1 < rev_map.size() && 
	rev_map.size() / 2 >= MIN_REV_MAP_SIZE)
    {
      rev_map.resize (rev_map.size() / 2);
    }

    // Get the current vertex index, and assign this index to the vertex
    // currently indexed (n - 1).
    unsigned int   index = index_map[v];

    if (index == n_vertices)
      return;
    
    Vertex_handle  last_v = rev_map[n_vertices];
    index_map[last_v] = index;
    rev_map[index] = last_v;

    // Clear the reverse mapping for the last vertex.
    rev_map[n_vertices] = Vertex_handle();
  }
  //@}

private:
  
  /*! Initialize the map for the given arrangement. */
  void _init ()
  {
    // Get the number of vertices and allocate the reverse map accordingly.
    n_vertices = static_cast<unsigned int>(this->arrangement()->number_of_vertices());
    
    if (n_vertices < MIN_REV_MAP_SIZE)
      rev_map.resize (MIN_REV_MAP_SIZE);
    else
      rev_map.resize (n_vertices);

    // Clear the current mapping.
    index_map.clear();

    // Create the initial mapping. 
    typename Arrangement_2::Vertex_iterator   vit;
    Vertex_handle                             vh;
    unsigned int                              index = 0;

    for (vit = this->arrangement()->vertices_begin();
	 vit != this->arrangement()->vertices_end(); ++vit, ++index)
    {
      // Map the current vertex to the current index.
      vh = vit;
      index_map[vh] = index;
      rev_map[index] = vh;
    }
  }  

};

/*!
 * Get the index property-map function. Provided so that boost is able to
 * access the Arr_vertex_index_map above.
 * \param index_map The index map.
 * \param v A vertex handle.
 * \return The vertex index.
 */
template<class Arrangement>
unsigned int get (const CGAL::Arr_vertex_index_map<Arrangement>& index_map,
		  typename Arrangement::Vertex_handle v) 
{ 
  return index_map[v];
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
