// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
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
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_K_NEIGHBOR_NEIGHBOR_H
#define CGAL_K_NEIGHBOR_NEIGHBOR_H

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_vertex_handle_3.h>
#include <list>

namespace CGAL {


/// Wrapper around Orthogonal_k_neighbor_search for Triangulation_3 Vertex_handles.
template <class Gt, class Vertex_handle>
class K_nearest_neighbor
{
public:
  typedef Gt  Geom_traits;

  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;
  typedef CGAL::Point_vertex_handle_3<Vertex_handle> Point_vertex_handle_3;
  typedef Search_traits_vertex_handle_3<Vertex_handle> Traits;
  typedef Euclidean_distance_vertex_handle_3<Vertex_handle> KDistance;
  typedef Orthogonal_k_neighbor_search<Traits,KDistance> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

private:
  Tree m_tree;

public:
  K_nearest_neighbor() {}

  /// @commentheading Precondition: 
  /// InputIterator value_type must be convertible to Point_vertex_handle_3.
  template <class InputIterator>
  K_nearest_neighbor(InputIterator first, InputIterator beyond)
  {
    m_tree = Tree(first, beyond);
  }

  /// Default copy constructor, operator =() and destructor are fine.

  /// @commentheading Precondition: 
  /// InputIterator value_type must be convertible to Point_vertex_handle_3.
  template <class InputIterator>
  void insert(InputIterator first, InputIterator beyond)
  {
    m_tree = Tree(first, beyond);
  }

  /// Empty KD-tree.
  void clear()
  {
    m_tree = Tree();
  }

  /// Search 'nb' nearest_neighbors of 'query' point.
  bool get_k_nearest_neighbors(const Point& query,
                               const unsigned int nb,
                               std::list<Vertex_handle>& nearest_neighbors)
  {
    Point_vertex_handle_3 point_wrapper(query.x(), query.y(), query.z(), NULL);
    Neighbor_search search(m_tree, point_wrapper, nb); // only nb nearest neighbors
    Search_iterator it = search.begin();
    for(unsigned int i=0;i<nb;i++,it++)
    {
      if(it == search.end())
        return false;
      nearest_neighbors.push_back((Vertex_handle)it->first);
    }
    return true;
  }

}; // K_nearest_neighbor


} //namespace CGAL

#endif // CGAL_K_NEIGHBOR_NEIGHBOR_H
