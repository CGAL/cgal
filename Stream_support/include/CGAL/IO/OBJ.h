// Copyright (c) 2015  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Lutz Kettner
//             Andreas Fabri
//             Maxime Gimeno

#ifndef CGAL_IO_OBJ_H
#define CGAL_IO_OBJ_H

#include <CGAL/IO/OBJ/File_writer_wavefront.h>
#include <CGAL/IO/Generic_writer.h>

#include <boost/range/value_type.hpp>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

namespace IO {
namespace internal {

template <typename PointRange, typename PolygonRange, typename VertexNormalOutputIterator>
bool read_OBJ(std::istream& is,
              PointRange& points,
              PolygonRange& faces,
              VertexNormalOutputIterator vn_out)
{
  typedef typename boost::range_value<PointRange>::type                               Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel                                 Kernel;
  typedef typename Kernel::Vector_3                                                   Normal;

  int mini(1),
      maxi(-INT_MAX);
  Point p;
  std::string line;

  while(getline(is, line))
  {
    if(line[0] == 'v' && line[1] == ' ')
    {
      std::istringstream iss(line.substr(1));
      iss >> p; // @fixme check successful reading
      if(!iss)
        return false;

      points.push_back(p);
    }
    else if(line[0] == 'f' && line[1] == ' ') // @fixme range checks
    {
      std::istringstream iss(line.substr(1));
      int i;
      faces.push_back(std::vector<std::size_t>());
      while(iss >> i)
      {
        if(i < 1)
        {
          faces.back().push_back(points.size()+i); // negative indices are relative references
          if(i < mini)
            mini = i;
        }
        else
        {
          faces.back().push_back(i-1);
          if(i-1 > maxi)
            maxi = i-1;
        }
        iss.ignore(256, ' ');
      }

      if(iss.fail())
        return false;
    }
    else if(line[0] == 'v' &&  line[1] == 'n' && line[2] == ' ')
    {
      std::istringstream iss(line.substr(1));
      double nx, ny, nz; // @fixme double?
      if(iss >> nx >> ny >> nz)
        *vn_out++ = Normal(nx, ny, nz); // @fixme check that every vertex has a normal?
      else
        return false;
    }
    else
    {
      //std::cerr << "ERROR : Cannnot read line beginning with " << line[0] << std::endl;
     continue;
    }
  }

  if(maxi > static_cast<int>(points.size()) || mini < -static_cast<int>(points.size()))
  {
    std::cerr << "a face index is invalid " << std::endl;
    return false;
  }

  return !is.fail();
}

} // namespace internal
} // namespace IO

//! \ingroup IOstreamFunctions
//!
/// reads the content of `is` into `points` and `faces`, using the `OBJ` format.
///
/// \tparam Points a `RandomAccessContainer` of `Point_3,
/// \tparam Faces a `RandomAccessContainer` of `RandomAccessContainer` of `std::size_t`
///
/// \see \ref IOStreamOBJ
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(std::istream& is,
              PointRange& points,
              PolygonRange& faces,
              const CGAL_BGL_NP_CLASS& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  return IO::internal::read_OBJ(is, points, faces,
                                choose_parameter(get_parameter(np, internal_np::vertex_normal_output_iterator),
                                                 CGAL::Emptyset_iterator()));
}

template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(const char* fname,
              PointRange& points,
              PolygonRange& polygons,
              const CGAL_BGL_NP_CLASS& np)
{
  std::ifstream in(fname);
  return read_OBJ(in, points, polygons, np);
}

template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(const std::string& fname, PointRange& points, PolygonRange& polygons, const CGAL_BGL_NP_CLASS& np)
{
  return read_OBJ(fname.c_str(), points, polygons, np);
}

template <typename PointRange, typename PolygonRange>
bool read_OBJ(std::istream& is, PointRange& points, PolygonRange& polygons)
{
  return read_OBJ(is, points, polygons, parameters::all_default());
}

template <typename PointRange, typename PolygonRange>
bool read_OBJ(const char* fname, PointRange& points, PolygonRange& polygons)
{
  return read_OBJ(fname, points, polygons, parameters::all_default());
}

template <typename PointRange, typename PolygonRange>
bool read_OBJ(const std::string& fname, PointRange& points, PolygonRange& polygons)
{
  return read_OBJ(fname, points, polygons, parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
 * \ingroup IOstreamFunctions
 *
 * writes the content of `points` and `polygons` in `os`, in the OBJ format.
 *
 * \see \ref IOStreamOBJ
 */
template <typename PointRange, typename PolygonRange>
bool write_OBJ(std::ostream& os,
               PointRange& points,
               PolygonRange& polygons)
{
  Generic_writer<std::ostream, File_writer_wavefront> writer(os);
  return writer(points, polygons);
}

template <typename PointRange, typename PolygonRange>
bool write_OBJ(const char* fname, PointRange& points, PolygonRange& polygons)
{
  std::ofstream out(fname);
  return write_OBJ(out, points, polygons);
}

template <typename PointRange, typename PolygonRange>
bool write_OBJ(const std::string& fname, PointRange& points, PolygonRange& polygons)
{
  return write_OBJ(fname.c_str(), points, polygons);
}

} // namespace CGAL

#endif // CGAL_IO_OBJ_H
