// Copyright (c) 2009-2013 GeometryFactory (France).
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
//
//
// Author(s)     : Laurent Rineau and Sebastien Loriot

#ifndef CGAL_POLYGON_SOUP_TO_POLYHEDRON_3_H
#define CGAL_POLYGON_SOUP_TO_POLYHEDRON_3_H

#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>

namespace CGAL{

namespace internal{
/**
  * Modifier to build a polyhedron from a soup of polygons.
  */
template <class HDS, class Point, class Polygon>
class Polygon_soup_to_polyhedron_3: public CGAL::Modifier_base<HDS>
{
  const std::vector<Point>& points;
  const std::vector<Polygon>& polygons;
  typedef typename HDS::Vertex::Point Point_3;
public:
  /**
   * The constructor for modifier object.
   * @param points points of the soup of polygons.
   * @param polygons each element in the vector describes a polygon using the index of the points in the vector.
   */
  Polygon_soup_to_polyhedron_3(
    const std::vector<Point>& points,
    const std::vector<Polygon>& polygons)
    : points(points), polygons(polygons)
  { }

  void operator()(HDS& out_hds)
  {
    Polyhedron_incremental_builder_3<HDS> builder(out_hds);

    builder.begin_surface(points.size(), polygons.size());

    for(std::size_t i = 0, end = points.size(); i < end; ++i)
      builder.add_vertex( Point_3(points[i][0], points[i][1], points[i][2]) );

    for(std::size_t i = 0, end = polygons.size(); i < end; ++i)
    {
      const Polygon& polygon = polygons[i];
      const std::size_t size = polygon.size();

      builder.begin_facet();
      for(std::size_t j = 0; j < size; ++j)
        builder.add_vertex_to_facet(polygon[j]);
      builder.end_facet();
    }
    builder.end_surface();
  }
};

} //namespace internal

/**
  * Append a soup of polygons in a Polyhedron
  */
template <class Polyhedron, class Point, class Polygon>
void polygon_soup_to_polyhedron_3(Polyhedron& P,
                                  const std::vector<Point>& points,
                                  const std::vector<Polygon>& polygons)
{
  internal::Polygon_soup_to_polyhedron_3< typename Polyhedron::HalfedgeDS,
                                          Point, Polygon > modifier(points, polygons);
  P.delegate(modifier);
}

} //end of namespace CGAL

#endif //CGAL_POLYGON_SOUP_TO_POLYHEDRON_3_H
