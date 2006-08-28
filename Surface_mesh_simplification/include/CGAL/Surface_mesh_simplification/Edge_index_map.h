// Copyright (c) 2006 Geometry Factory (France).
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Surface_mesh_simplification/Polyhedron_edge_cached_pointer_map.h $
// $Id: Polyhedron_edge_cached_pointer_map.h 32048 2006-06-23 13:59:36Z lsaboret $
// 
//
// Author(s): Fernando Cacciola <fernando.cacciola@gmail.com>

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_STORED_INDEX_MAP_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_STORED_INDEX_MAP_H

namespace boost
{

template<class Graph>
class Edge_stored_index_map : public put_get_helper< typename Graph_::size_type, Edge_stored_index_map<Graph_> >
{
private:

  typedef Graph_ Graph ;

public:

  typedef typename G::size_type size_type ;
  
  typedef readable_property_map_tag                           category;
  typedef size_type                                           value_type;
  typedef size_type                                           reference;
  typedef typename graph_traits<Graph const>::edge_descriptor key_type;

  reference operator[](key_type const& e) const { return e->id(); }
};

} // namespace boost

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_STORED_INDEX_MAP_H
