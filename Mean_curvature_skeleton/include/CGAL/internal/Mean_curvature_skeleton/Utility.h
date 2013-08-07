// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
//
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_MCFSKEL_UTILITY_H
#define CGAL_MCFSKEL_UTILITY_H

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <boost/graph/graph_traits.hpp>

namespace CGAL {
namespace internal {

template<class Polyhedron>
typename boost::graph_traits<Polyhedron>::edge_descriptor
mesh_split(Polyhedron& polyhedron, 
           typename boost::graph_traits<Polyhedron>::edge_descriptor ei,
           typename Polyhedron::Traits::Point_3 pn)
{
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor            edge_descriptor;

  edge_descriptor en = polyhedron.split_edge(ei);
  en->vertex()->point() = pn;
  polyhedron.split_facet(en, ei->next());

  en->id() = -1;
  en->opposite()->id() = -1;
  ei->id() = -1;
  ei->opposite()->id() = -1;
  en->next()->id() = -1;
  en->next()->opposite()->id() = -1;
  en->next()->next()->id() = -1;
  ei->next()->id() = -1;
  edge_descriptor ej = en->opposite();
  if (!(ej->is_border()))
  {
    polyhedron.split_facet(ei->opposite(), ej->next());
    ej->next()->id() = -1;
    edge_descriptor ei_op_next = ei->opposite()->next();
    ei_op_next->id() = -1;
    ei_op_next->opposite()->id() = -1;
    ei_op_next->next()->id() = -1;
  }

  return en;
}

template<class Vertex, class Kernel>
double get_triangle_area(Vertex v1,
                         Vertex v2,
                         Vertex v3)
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  Point p1 = v1->point();
  Point p2 = v2->point();
  Point p3 = v3->point();
  Vector v12(p1, p2);
  Vector v13(p1, p3);
  return sqrtf(cross_product(v12, v13).squared_length()) * 0.5;
}

template<class Polyhedron>
double get_surface_area(Polyhedron& polyhedron)
{
  typedef typename Polyhedron::Traits                                 Kernel;
  typedef typename Polyhedron::Facet_iterator                         Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator       Halfedge_facet_circulator;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor	vertex_descriptor;

  double total_area = 0;
  for (Facet_iterator i = polyhedron.facets_begin(); i != polyhedron.facets_end(); ++i)
  {
    Halfedge_facet_circulator j = i->facet_begin();
    vertex_descriptor v1 = j->vertex();
    ++j;
    vertex_descriptor v2 = j->vertex();
    ++j;
    vertex_descriptor v3 = j->vertex();
    total_area += internal::get_triangle_area<vertex_descriptor, Kernel>(v1, v2, v3);
  }
  return total_area;
}

} //namespace internal
} //namespace CGAL

#endif //CGAL_MCFSKEL_UTILITY_H
