// Copyright (c) 2011 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot


#ifndef CGAL_TRIANGLE_ACCESSOR_WITH_PPMAP_3_H
#define CGAL_TRIANGLE_ACCESSOR_WITH_PPMAP_3_H

#include <CGAL/license/Polygon_mesh_processing.h>


#include <CGAL/Kernel_traits.h>
#include <CGAL/property_map.h>

namespace CGAL{

template < class Polyhedron, class PolyhedronPointPMap>
struct Triangle_accessor_with_ppmap_3
{
  typedef PolyhedronPointPMap Ppmap;
  typedef typename boost::property_traits< Ppmap >::value_type          Point_3;
  typedef typename CGAL::Kernel_traits<Point_3>::Kernel                  Kernel;
  typedef typename Kernel::Triangle_3                                Triangle_3;
  typedef typename Polyhedron::Facet_const_iterator           Triangle_iterator;
  typedef typename Polyhedron::Facet_const_handle               Triangle_handle;

  PolyhedronPointPMap ppmap;

  Triangle_accessor_with_ppmap_3(){}

  Triangle_accessor_with_ppmap_3(PolyhedronPointPMap ppmap):ppmap(ppmap) {}

  Triangle_iterator triangles_begin(const Polyhedron& p) const
  {
    return p.facets_begin();
  }

  Triangle_iterator triangles_end(const Polyhedron& p) const
  {
    return p.facets_end();
  }

  Triangle_3 triangle(const Triangle_handle& handle) const
  {
    typedef typename Kernel::Point_3 Point;
    const Point& a = get(ppmap, handle->halfedge()->vertex());
    const Point& b = get(ppmap, handle->halfedge()->next()->vertex());
    const Point& c = get(ppmap, handle->halfedge()->next()->next()->vertex());
    return Triangle_3(a,b,c);
  }
};

} // end of namespace CGAL

#endif // CGAL_TRIANGLE_ACCESSOR_WITH_PPMAP_3_H
