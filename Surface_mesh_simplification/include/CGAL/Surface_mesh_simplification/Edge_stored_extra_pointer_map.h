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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Surface_mesh_simplification/Edge_stored_extra_pointer_map.h $
// $Id: Edge_stored_extra_pointer_map.h 32048 2006-06-23 13:59:36Z lsaboret $
// 
//
// Author(s): Fernando Cacciola <fernando.cacciola@gmail.com>

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_STORED_EXTRA_POINTER_MAP_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_STORED_EXTRA_POINTER_MAP_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

CGAL_BEGIN_NAMESPACE

template<class Graph_>
class Edge_stored_extra_pointer_map : public boost::put_get_helper<void*&, Edge_stored_extra_pointer_map<Graph_> >
{
private:

  typedef Graph_ Graph ;
  
public:

  typedef boost::lvalue_property_map_tag                       category;
  typedef void*                                                value_type;
  typedef void*&                                               reference;
  typedef typename boost::graph_traits<Graph>::edge_descriptor key_type;

  reference operator[](key_type const& e) const { return e->extra_pointer() ; }
};

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_STORED_EXTRA_POINTER_MAP_H
