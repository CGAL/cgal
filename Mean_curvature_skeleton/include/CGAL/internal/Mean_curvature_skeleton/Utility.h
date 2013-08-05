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
#ifndef CGAL_MCFSKELETON_UTILITY_H
#define CGAL_MCFSKELETON_UTILITY_H

namespace CGAL {
namespace internal {

template<class Polyhedron, class edge_descriptor, class Point>
edge_descriptor mesh_split(Polyhedron& polyhedron, edge_descriptor ei, Point pn)
{
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

} //namespace internal
} //namespace CGAL

#endif //CGAL_MCFSKELETON_UTILITY_H
