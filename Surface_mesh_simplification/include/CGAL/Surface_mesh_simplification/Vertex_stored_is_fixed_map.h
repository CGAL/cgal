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
// $URL$
// $Id$
// 
//
// Author(s): Fernando Cacciola <fernando.cacciola@gmail.com>

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_STORED_IS_VERTEX_FIXED_MAP_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_STORED_IS_VERTEX_FIXED_MAP_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

CGAL_BEGIN_NAMESPACE

template<class Graph_>
class Vertex_stored_is_fixed_map : public boost::put_get_helper<bool, Vertex_stored_is_fixed_map<Graph_> >
{
private:

  typedef Graph_ Graph ;
  
public:

  typedef boost::readable_property_map_tag                             category;
  typedef bool                                                         value_type;
  typedef bool                                                         reference;
  typedef typename boost::graph_traits<Graph const>::vertex_descriptor key_type;

  reference operator[](key_type const& v) const { return v->is_fixed() ; }
};

        
CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_STORED_IS_VERTEX_FIXED_MAP_H
