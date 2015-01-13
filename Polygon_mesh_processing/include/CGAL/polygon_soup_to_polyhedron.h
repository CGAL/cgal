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
// Author(s)     : Laurent Rineau and Ilker O. Yaz

#ifndef CGAL_POLYGON_SOUP_TO_POLYHEDRON
#define CGAL_POLYGON_SOUP_TO_POLYHEDRON

#include <CGAL/Modifier_base.h>
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace CGAL
{
/**
* \ingroup PkgPolygonMeshProcessing
* Polyhedron modifier to build a polyhedron from a soup of polygons.
* \todo convert to a function that is creating a model of `MutableFaceGraph`
*/
template <class HDS, class Point_3>
class Polygon_soup_to_polyhedron_3 : public CGAL::Modifier_base<HDS>
{
  typedef std::vector<std::size_t> Polygon_3;

  const std::vector<Point_3>& points;
  const std::vector<std::vector<std::size_t> >& polygons;
public:
  /**
  * The constructor for modifier object.
  * @param points points of the soup of polygons.
  * @param polygons each element in the vector describes a polygon using the index of the points in the vector.
  */
  Polygon_soup_to_polyhedron_3(const std::vector<Point_3>& points,
    const std::vector<std::vector<std::size_t> >& polygons)
    : points(points), polygons(polygons)
  { }

  void operator()(HDS& out_hds)
  {
    Polyhedron_incremental_builder_3<HDS> builder(out_hds);

    builder.begin_surface(points.size(), polygons.size());

    for (std::size_t i = 0, end = points.size(); i < end; ++i)
    {
      builder.add_vertex(points[i]);
    }

    for (std::size_t i = 0, end = polygons.size(); i < end; ++i)
    {
      const Polygon_3& polygon = polygons[i];
      const std::size_t size = polygon.size();

      builder.begin_facet();
      for (std::size_t j = 0; j < size; ++j) {
        builder.add_vertex_to_facet(polygon[j]);
      }
      builder.end_facet();
    }
    builder.end_surface();
  }
};

}// end namespace CGAL

#endif
