// Copyright (c) 2015 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau and Sebastien Loriot

#ifndef CGAL_IO_OFF_H
#define CGAL_IO_OFF_H

#include <CGAL/IO/reader_helpers.h>
#include <CGAL/IO/OFF/File_scanner_OFF.h>
#include <CGAL/IO/OFF/File_writer_OFF.h>

#include <vector>
#include <iostream>
#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>

namespace CGAL {

/*!
 * \ingroup IOstreamFunctions
 *
 * reads the content of `in` into `points` and `polygons`, in the COFF format.
 *
 * \see \ref IOStreamOFF
 */
template <class Point_3, class Polygon_3, class Color_rgb >
bool read_OFF(std::istream& in,
              std::vector< Point_3 >& points,
              std::vector< Polygon_3 >& polygons,
              std::vector<Color_rgb>& fcolors,
              std::vector<Color_rgb>& vcolors,
              bool /* verbose */ = false)
{
  CGAL::File_scanner_OFF scanner(in);

  points.resize(scanner.size_of_vertices());
  polygons.resize(scanner.size_of_facets());

  if(scanner.has_colors())
    vcolors.resize(scanner.size_of_vertices());

  for(std::size_t i = 0; i < scanner.size_of_vertices(); ++i)
  {
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    CGAL_assertion(w!=0);
    IO::internal::fill_point( x/w, y/w, z/w, points[i] );
    if(scanner.has_colors())
    {
      unsigned char r=0, g=0, b=0;
      scanner.scan_color( r, g, b);
      vcolors[i] = Color_rgb(r,g,b);
    }
    else
    {
      scanner.skip_to_next_vertex(i);
    }

    if(!in)
      return false;
  }

  bool has_fcolors = false;
  for(std::size_t i = 0; i < scanner.size_of_facets(); ++i)
  {
    std::size_t no;
    scanner.scan_facet( no, i);

    if(!in)
      return false;

    IO::internal::resize(polygons[i], no);
    for(std::size_t j = 0; j < no; ++j)
    {
      std::size_t id;
      scanner.scan_facet_vertex_index(id, i);
      if(id < scanner.size_of_vertices())
        polygons[i][j] = id;
      else
        return false;
    }

    if(i == 0)
    {
      std::string col;
      std::getline(in, col);
      std::istringstream iss(col);
      char ci =' ';

      if(iss >> ci)
      {
        has_fcolors = true;
        fcolors.resize(scanner.size_of_facets());
        std::istringstream iss2(col);
        fcolors[i] = scanner.get_color_from_line(iss2);
      }
    }
    else if(has_fcolors)
    {
      unsigned char r=0, g=0, b=0;
      scanner.scan_color(r,g,b);
      fcolors[i] = Color_rgb(r,g,b);
    }
  }

  return in.good();
}

/*!
 * \ingroup IOstreamFunctions
 *
 * reads the content of `in` into `points` and `polygons`, in the OFF format.
 *
 * \see \ref IOStreamOFF
 */
template <class Point_3, class Polygon_3>
bool read_OFF(std::istream& in,
              std::vector< Point_3 >& points,
              std::vector< Polygon_3 >& polygons,
              bool /* verbose */ = false)
{
  CGAL::File_scanner_OFF scanner(in);

  points.resize(scanner.size_of_vertices());
  polygons.resize(scanner.size_of_facets());
  for(std::size_t i = 0; i < scanner.size_of_vertices(); ++i)
  {
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    CGAL_assertion(w!=0);
    IO::internal::fill_point( x/w, y/w, z/w, points[i] );
    scanner.skip_to_next_vertex( i);
  }

  if(!in)
    return false;

  for(std::size_t i = 0; i < scanner.size_of_facets(); ++i)
  {
    std::size_t no;

    scanner.scan_facet( no, i);
    IO::internal::resize(polygons[i], no);
    for(std::size_t j = 0; j < no; ++j)
    {
      std::size_t id;
      scanner.scan_facet_vertex_index(id, i);
      if(id < scanner.size_of_vertices())
        polygons[i][j] = id;
      else
        return false;
    }
  }

  return in.good();
}

/*!
 * \ingroup IOstreamFunctions
 *
 * writes the content of `points` and `polygons` in `out`, in the OFF format.
 *
 * \see \ref IOStreamOFF
 */
template <class Point_3, class Polygon_3>
bool write_OFF(std::ostream& out,
               std::vector< Point_3 >& points,
               std::vector< Polygon_3 >& polygons)
{
  CGAL::File_writer_OFF writer;
  writer.write_header(out, points.size(), 0, polygons.size());

  for(std::size_t i = 0, end = points.size(); i < end; ++i)
  {
    const Point_3& p = points[i];
    writer.write_vertex( p.x(), p.y(), p.z() );
  }

  writer.write_facet_header();
  for(std::size_t i=0, end=polygons.size(); i<end; ++i)
  {
    Polygon_3& polygon = polygons[i];
    const std::size_t size = polygon.size();

    writer.write_facet_begin(size);
    for(std::size_t j=0; j<size; ++j)
      writer.write_facet_vertex_index(polygon[j]);
    writer.write_facet_end();
  }
  writer.write_footer();

  return out.good();
}

} // namespace CGAL

#endif // CGAL_IO_OFF_H
