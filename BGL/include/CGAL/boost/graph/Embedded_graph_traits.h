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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Polyhedron_extended_BGL.h $
// $Id: Polyhedron_extended_BGL.h 32155 2006-06-29 17:08:41Z fcacciola $
// 
//
// Author(s) : Fernando Caccciola <fernando.cacciola@gmail.com>

#ifndef CGAL_BOOST_GRAPH_GEOMETRIC_GRAPH_TRAITS_H
#define CGAL_BOOST_GRAPH_GEOMETRIC_GRAPH_TRAITS_H

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map.hpp>

#include <CGAL/basic.h>
#include <CGAL/boost/graph/BGL_properties.h>

CGAL_BEGIN_NAMESPACE

template<class Graph> struct Geometric_graph_traits ;

template<class Graph>
inline
typename Geometric_graph_traits<Graph>::Point const& 
get_point ( typename boost::graph_traits<Graph const>::vertex_descriptor const& v, Graph const& g )
{
  return boost::get(vertex_point_t(),g,v) ;
}


template<class Graph>
inline
void set_point ( typename boost::graph_traits<Graph>::vertex_descriptor const& v
               , Graph&                                                        g
               , typename Geometric_graph_traits<Graph>::Point const&          p
               )
{
  boost::put(vertex_point_t(),g,v,p) ;
}

CGAL_END_NAMESPACE

#endif // CGAL_BOOST_GRAPH_GEOMETRIC_GRAPH_TRAITS_H
