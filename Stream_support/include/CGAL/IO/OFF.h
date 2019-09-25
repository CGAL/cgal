// Copyright (c) 2015 GeometryFactory
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Laurent Rineau and Sebastien Loriot

#ifndef CGAL_IO_OFF_H
#define CGAL_IO_OFF_H

#include <CGAL/IO/OFF/OFF_reader.h>
#include <CGAL/IO/OFF/File_scanner_OFF.h>
#include <CGAL/IO/reader_helpers.h>
#include <CGAL/IO/OFF/File_writer_OFF.h>

#include <vector>
#include <iostream>
#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>

namespace CGAL {
template <class Point_3, class Polygon_3>
bool
read_OFF( std::istream& in,
          std::vector< Point_3 >& points,
          std::vector< Polygon_3 >& polygons,
          bool /* verbose */ = false)
{
  return OFF_internal::read_OFF(in, points, polygons);
}

template <class Point_3, class Polygon_3, class Color_rgb >
bool
read_OFF( std::istream& in,
          std::vector< Point_3 >& points,
          std::vector< Polygon_3 >& polygons,
          std::vector<Color_rgb>& fcolors,
          std::vector<Color_rgb>& vcolors,
          bool /* verbose */ = false)
{
  return OFF_internal::read_OFF(in, points, polygons, fcolors, vcolors);
}

template <class Point_3, class Polygon_3>
bool
write_OFF(std::ostream& out,
          std::vector< Point_3 >& points,
          std::vector< Polygon_3 >& polygons)
{
  CGAL::File_writer_OFF writer;
  writer.write_header(out,
                      points.size(),
                      0,
                      polygons.size());
  for(std::size_t i = 0, end = points.size();
      i < end; ++i)
  {
    const Point_3& p = points[i];
    writer.write_vertex( p.x(), p.y(), p.z() );
  }
  writer.write_facet_header();
  for(std::size_t i = 0, end = polygons.size();
      i < end; ++i)
  {
    Polygon_3& polygon = polygons[i];
    const std::size_t size = polygon.size();
    writer.write_facet_begin(size);
    for(std::size_t j = 0; j < size; ++j) {
      writer.write_facet_vertex_index(polygon[j]);
    }
    writer.write_facet_end();
  }
  writer.write_footer();

  return (bool) out;
}

} // namespace CGAL

#endif // CGAL_IO_OFF_H
