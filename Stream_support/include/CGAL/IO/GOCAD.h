// Copyright (c) 2020 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//                 Mael Rouxel-Labbé
//                 Maxime Gimeno

#ifndef CGAL_IO_GOCAD_H
#define CGAL_IO_GOCAD_H

#include <CGAL/IO/reader_helpers.h>
#include <CGAL/IO/Generic_writer.h>

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/use.h>

#include <boost/range/value_type.hpp>

#include <fstream>
#include <iostream>
#include <vector>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read
template <typename PointRange, typename PolygonRange, typename NamedParameters>
bool read_GOCAD(std::istream& input,
                std::pair<std::string, std::string>& name_and_color,
                PointRange& points,
                PolygonRange& polygons,
                const NamedParameters&)
{
  typedef typename boost::range_value<PointRange>::type     Point;
  typedef typename boost::range_value<PolygonRange>::type   CGAL_Polygon;
  int offset = 0;
  char c;
  std::string s, tface("TFACE");
  int i,j,k;
  Point p;
  bool vertices_read = false;

  while(input >> s)
  {
    if(s == tface)
      break;

    std::string::size_type idx;
    if((idx = s.find("name")) != std::string::npos)
    {
      std::istringstream str(s.substr(idx + 5));
      str >> name_and_color.first;
    }

    if((idx = s.find("color")) != std::string::npos)
    {
      std::istringstream str(s.substr(idx + 6));
      str >> name_and_color.second;
    }
  }
  std::getline(input, s);

  while(input.get(c))
  {
    if((c == 'V') || (c == 'P'))
    {
      input >> s >> i >> p; // @fixme check for failure
      if(!vertices_read)
      {
        vertices_read = true;
        offset -= i; // Some files start with index 0 others with 1
      }

      points.push_back(p);
    }
    else if(vertices_read && (c == 'T'))
    {
      input >> c >> c >> c >>  i >> j >> k;
      CGAL_Polygon new_face(3);
      new_face[0] = offset+i;
      new_face[1] = offset+j;
      new_face[2] = offset+k;
      polygons.push_back(new_face);
    }
    else if(c == 'E')
    {
      break;
    }

    std::getline(input, s);
  }

  return !input.fail();
}
/*!
 * \ingroup IOstreamFunctions
 *
 * reads the content of `is` into `points` and `polygons`, in the GOCAD format.
 *
 * \see \ref IOStreamOFF
 */
template <typename PointRange, typename PolygonRange, typename NamedParameters>
bool read_GOCAD(std::istream& is,
                PointRange& points,
                PolygonRange& polygons,
                const NamedParameters&np)
{
  std::pair<std::string, std::string> dummy;
  return read_GOCAD(is, dummy, points, polygons, np);
}

/*!
 * \ingroup IOstreamFunctions
 *
 * reads the content of the file `fname` into `points` and `polygons`, in the GOCAD format.
 *
 * \see \ref IOStreamOFF
 */
template <typename PointRange, typename PolygonRange, typename NamedParameters>
bool read_GOCAD(const char* fname,
                PointRange& points,
                PolygonRange& polygons,
                const NamedParameters& np)
{
  std::ifstream in(fname);
  std::pair<std::string, std::string> dummy;
  return read_GOCAD(in, dummy, points, polygons, np);
}



template <typename PointRange, typename PolygonRange>
bool read_GOCAD(std::istream& is, PointRange& points, PolygonRange& polygons)
{
  std::pair<std::string, std::string> dummy;
  return read_GOCAD(is, dummy, points, polygons, parameters::all_default());
}

template <typename PointRange, typename PolygonRange>
bool read_GOCAD(const char* fname, PointRange& points, PolygonRange& polygons)
{
  return read_GOCAD(fname, points, polygons, parameters::all_default());
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write
namespace IO{
namespace internal{

template <typename PointRange,
          typename PolygonRange,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(std::ostream& os,
                 const char* fname,
                 const PointRange& points,
                 const PolygonRange& polygons,
                 const CGAL_BGL_NP_CLASS&)
{
  typedef typename boost::range_value<PointRange>::type     Point;
  typedef typename boost::range_value<PolygonRange>::type   CGAL_polygon;
  if(!os.good())
    return false;

  os << "GOCAD TSurf 1\n"
        "HEADER {\n"
        "name:";
  os << fname << std::endl;
  os << "*border:on\n"
        "*border*bstone:on\n"
        "}\n"
        "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n"
        "NAME Default\n"
        "AXIS_NAME \"X\" \"Y\" \"Z\"\n"
        "AXIS_UNIT \"m\" \"m\" \"m\"\n"
        "ZPOSITIVE Elevation\n"
        "END_ORIGINAL_COORDINATE_SYSTEM\n"
        "TFACE\n";

  std::size_t i = 0;
  for(const Point& p : points)
  {
    os << "VRTX " << i << " " << p << "\n";
  }

  for(const CGAL_polygon& poly : polygons)
  {
    os << "TRGL";
    for(const auto& id : poly)
      os << " " << points[id];
    os<< "\n";
  }

  os << "END" << std::endl;

  return os.good();
}
} //end internal
}//end IO
/*!
  \ingroup IOstreamFunctions

 * writes the content of `points` and `polygons` in `out`, in the TS format.

  \see \ref IOStreamGocad
*/
template <typename PointRange,
          typename PolygonRange,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(std::ostream& os,
                 const PointRange& points,
                 const PolygonRange& polygons,
                 const CGAL_BGL_NP_CLASS&np)
{
  return IO::internal::write_GOCAD(os, "anonymous", points, polygons, np);

}

template <typename PointRange,
          typename PolygonRange>
bool write_GOCAD(std::ostream& os,
                 const PointRange& points,
                 const PolygonRange& polygons)
{
  return IO::internal::write_GOCAD(os, "anonymous", points, polygons, parameters::all_default());

}

template <typename PointRange,
          typename PolygonRange,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(const char* fname,
                 const PointRange& points,
                 const PolygonRange& polygons,
                 const CGAL_BGL_NP_CLASS&np)
{
  std::ofstream os(fname);
  return IO::internal::write_GOCAD(os, fname, points, polygons, np);
}

template <typename PointRange,
          typename PolygonRange>
bool write_GOCAD(const char* fname,
                 const PointRange& points,
                 const PolygonRange& polygons)
{
  return IO::internal::write_GOCAD(fname, points, polygons, parameters::all_default());
}

template <typename PointRange,
          typename PolygonRange,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(const std::string& fname,
                 const PointRange& points,
                 const PolygonRange& polygons,
                 const CGAL_BGL_NP_CLASS&np)
{
  std::ofstream os(fname.c_str());
  return IO::internal::write_GOCAD(os, fname, points, polygons, np);
}

template <typename PointRange,
          typename PolygonRange>
bool write_GOCAD(const std::string& fname,
                 const PointRange& points,
                 const PolygonRange& polygons)
{
  return write_GOCAD(fname, points, polygons, parameters::all_default());
}


} // namespace CGAL

#endif // CGAL_IO_GOCAD_H
