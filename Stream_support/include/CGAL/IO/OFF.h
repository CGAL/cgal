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

#include <CGAL/IO/OFF/File_scanner_OFF.h>
#include <CGAL/IO/OFF/File_writer_OFF.h>
#include <CGAL/IO/reader_helpers.h>
#include <CGAL/IO/Generic_writer.h>

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/use.h>

#include <fstream>
#include <iostream>
#include <vector>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

namespace IO {
namespace internal {

template <typename Point_3, typename Polygon_3,
          typename VertexNormalOutputIterator,
          typename VertexColorOutputIterator,
          typename VertexTextureOutputIterator,
          typename FaceColorOutputIterator>
bool read_OFF(std::istream& in,
              std::vector<Point_3>& points,
              std::vector<Polygon_3>& polygons,
              VertexNormalOutputIterator vn_out,
              VertexColorOutputIterator vc_out,
              VertexTextureOutputIterator vt_out,
              FaceColorOutputIterator fc_out)

{
  typedef typename CGAL::Kernel_traits<Point_3>::Kernel                               Kernel;
  typedef typename Kernel::Point_2                                                    Texture;
  typedef typename Kernel::Vector_3                                                   Normal;
  typedef CGAL::Color                                                                 Color;

  CGAL::File_scanner_OFF scanner(in);

  points.resize(scanner.size_of_vertices());
  polygons.resize(scanner.size_of_facets());

  for(std::size_t i=0; i<scanner.size_of_vertices(); ++i)
  {
    double x, y, z, w;
    scanner.scan_vertex(x, y, z, w);
    CGAL_assertion(w != 0);
    IO::internal::fill_point(x, y, z, w, points[i]);

    if(scanner.has_normals())
    {
      double nx, ny, nz, nw;
      scanner.scan_normal(nx, ny, nz, nw);
      CGAL_assertion(nw != 0);
      *vn_out++ = Normal(nx, ny, nz, nw);
    }

    if(scanner.has_colors())
    {
      unsigned char r=0, g=0, b=0;
      scanner.scan_color(r, g, b);
      *vc_out++ = Color(r,g,b);
    }

    // @fixme
//    if(scanner.has_textures())
//    {
//      double nx, ny, nw;
//      scanner.scan_texture(nx, ny, nw);
//      CGAL_assertion(nw != 0);
//      *vt_out++ = Texture(nx, ny, nw);
//    }

    scanner.skip_to_next_vertex(i);

    if(!in.good())
      return false;
  }

  bool has_fcolors = false;
  for(std::size_t i=0; i<scanner.size_of_facets(); ++i)
  {
    std::size_t no;
    scanner.scan_facet(no, i);

    if(!in.good())
      return false;

    IO::internal::resize(polygons[i], no);
    for(std::size_t j=0; j<no; ++j)
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
        std::istringstream iss2(col);
        *fc_out++ = scanner.get_color_from_line(iss2);
      }
    }
    else if(has_fcolors)
    {
      unsigned char r=0, g=0, b=0;
      scanner.scan_color(r,g,b);
      *fc_out++ = Color(r,g,b);
    }
  }

  return !in.fail();
}

} // namespace internal
} // namespace IO

/*!
 * \ingroup IOstreamFunctions
 *
 * reads the content of `in` into `points` and `polygons`, in the OFF format.
 *
 * @fixme re-add support for STOFF
 *
 * \see \ref IOStreamOFF
 */
template <typename Point_3, typename Polygon_3, typename NamedParameters>
bool read_OFF(std::istream& in,
              std::vector<Point_3>& points, // @fixme doesn't need a vector specifically
              std::vector<Polygon_3>& polygons,
              const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  return IO::internal::read_OFF(in, points, polygons,
                                choose_parameter(get_parameter(np, internal_np::vertex_normal_output_iterator),
                                                 CGAL::Emptyset_iterator()),
                                choose_parameter(get_parameter(np, internal_np::vertex_color_output_iterator),
                                                 CGAL::Emptyset_iterator()),
                                choose_parameter(get_parameter(np, internal_np::vertex_texture_output_iterator),
                                                 CGAL::Emptyset_iterator()),
                                choose_parameter(get_parameter(np, internal_np::face_color_output_iterator),
                                                 CGAL::Emptyset_iterator()));
}

template <typename Point_3, typename Polygon_3, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const char* fname,
              std::vector<Point_3>& points,
              std::vector<Polygon_3>& polygons,
              const CGAL_BGL_NP_CLASS& np)
{
  std::ifstream in(fname);
  return read_OFF(in, points, polygons, np);
}

template <typename Point_3, typename Polygon_3, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const std::string& fname,
              std::vector<Point_3>& points,
              std::vector<Polygon_3>& polygons,
              const CGAL_BGL_NP_CLASS& np)
{
  return read_OFF(fname.c_str(), points, polygons, np);
}

template <typename Point_3, typename Polygon_3>
bool read_OFF(std::istream& is,
              std::vector<Point_3>& points,
              std::vector<Polygon_3>& polygons)
{
  return read_OFF(is, points, polygons, parameters::all_default());
}

template <typename Point_3, typename Polygon_3>
bool read_OFF(const char* fname,
              std::vector<Point_3>& points,
              std::vector<Polygon_3>& polygons)
{
  return read_OFF(fname, points, polygons, parameters::all_default());
}

template <typename Point_3, typename Polygon_3>
bool read_OFF(const std::string& fname,
              std::vector<Point_3>& points,
              std::vector<Polygon_3>& polygons)
{
  return read_OFF(fname, points, polygons, parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
 * \ingroup IOstreamFunctions
 *
 * writes the content of `points` and `polygons` in `out`, in the OFF format.
 *
 * \see \ref IOStreamOFF
 */
template <typename Point_3, typename Polygon_3, typename NamedParameters>
bool write_OFF(std::ostream& os,
               std::vector<Point_3>& points,
               std::vector<Polygon_3>& polygons,
               const NamedParameters& np)
{
  Generic_writer<std::ostream, File_writer_OFF> writer(os);
  return writer(points, polygons, np);
}

template <class Point_3, class Polygon_3>
bool write_OFF(std::ostream& os,
               std::vector<Point_3>& points,
               std::vector<Polygon_3>& polygons)
{
  return write_OFF(os, points, polygons, parameters::all_default());
}

} // namespace CGAL

#endif // CGAL_IO_OFF_H
