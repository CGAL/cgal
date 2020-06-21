// Copyright (c) 2016 GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_IO_XYZ_H
#define CGAL_POINT_SET_IO_XYZ_H

#include <CGAL/license/Point_set_3.h>

#include <CGAL/Point_set_3.h>

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <fstream>
#include <string>

namespace CGAL {

template <typename Point, typename Vector>
class Point_set_3;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

/*!
  \ingroup PkgPointSet3IO

  Reads the content of an intput stream in the XYZ format into a point set.

  \tparam Point a `CGAL::Point_3`
  \tparam Vector a `CGAL::Vector_3`

  \param stream the input stream
  \param point_set the point set

  \return `true` if the reading was successful, `false` otherwise.

  \see \ref IOStreamXYZ
 */
template <typename Point, typename Vector>
bool read_XYZ(std::istream& stream,
              CGAL::Point_set_3<Point, Vector>& point_set)
{
  point_set.add_normal_map();

  bool out = CGAL::read_xyz_points(stream,
                                   point_set.index_back_inserter(),
                                   CGAL::parameters::point_map(point_set.point_push_map())
                                                    .normal_map(point_set.normal_push_map()));

  bool has_normals = false;
  for(typename CGAL::Point_set_3<Point, Vector>::const_iterator it=point_set.begin(); it!=point_set.end(); ++it)
  {
    if(point_set.normal(*it) != CGAL::NULL_VECTOR)
    {
      has_normals = true;
      break;
    }
  }

  if(!has_normals)
    point_set.remove_normal_map();

  return out;
}

/*!
  \ingroup PkgPointSet3IO

  Reads the content of an input XYZ file in a point set.

  \tparam Point a `CGAL::Point_3`
  \tparam Vector a `CGAL::Vector_3`

  \param fname the path to the input file
  \param point_set the point set

  \return `true` if the reading was successful, `false` otherwise.

  \see \ref IOStreamXYZ
*/
template <typename Point, typename Vector>
bool read_XYZ(const char* fname, CGAL::Point_set_3<Point, Vector>& point_set)
{
  std::ifstream is(fname);
  return read_XYZ(is, point_set);
}

template <typename Point, typename Vector>
bool read_XYZ(const std::string& fname, CGAL::Point_set_3<Point, Vector>& point_set)
{
  return read_XYZ(fname.c_str(), point_set);
}

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgPointSet3IODeprecated

  \deprecated This function is deprecated since \cgal 5.2, `CGAL::read_XYZ()` should be used instead.
 */
template <typename Point, typename Vector>
CGAL_DEPRECATED bool read_xyz_point_set(std::istream& stream, ///< input stream.
                                        CGAL::Point_set_3<Point, Vector>& point_set) ///< point set
{
  return read_XYZ(stream, point_set);
}

#endif // CGAL_NO_DEPRECATED_CODE

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
  \ingroup PkgPointSet3IO

  Writes the content of a point set into an output stream in the XYZ format.

  \tparam Point a `CGAL::Point_3`
  \tparam Vector a `CGAL::Vector_3`

  \param stream the output stream
  \param point_set the point set.

  \return `true` if the writing was successful, `false` otherwise.

  \see \ref IOStreamXYZ
 */
template <typename Point, typename Vector>
bool write_XYZ(std::ostream& stream,
               const CGAL::Point_set_3<Point, Vector>& point_set)
{
  if(point_set.has_normal_map())
    return CGAL::write_XYZ(stream, point_set,
                           CGAL::parameters::point_map(point_set.point_map())
                                            .normal_map(point_set.normal_map()));

  return CGAL::write_XYZ(stream, point_set, CGAL::parameters::point_map(point_set.point_map()));
}

/*!
  \ingroup PkgPointSet3IO

  Writes the content of a point set into an output file in the XYZ format.

  \tparam Point a `CGAL::Point_3`
  \tparam Vector a `CGAL::Vector_3`

  \param fname the path to the output file
  \param point_set the point set.

  \return `true` if the writing was successful, `false` otherwise.

  \see \ref IOStreamXYZ
 */
template <typename Point, typename Vector>
bool write_XYZ(const char* fname, const CGAL::Point_set_3<Point, Vector>& point_set)
{
  std::ofstream os(fname);
  return write_XYZ(os, point_set);
}

template <typename Point, typename Vector>
bool write_XYZ(const std::string& fname, const CGAL::Point_set_3<Point, Vector>& point_set)
{
  return write_XYZ(fname.c_str(), point_set);
}

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgPointSet3IODeprecated

  \deprecated This function is deprecated since \cgal 5.2, `CGAL::write_XYZ()` should be used instead.
 */
template <typename Point, typename Vector>
CGAL_DEPRECATED bool write_xyz_point_set(std::ostream& stream,
                                         const CGAL::Point_set_3<Point, Vector>& point_set)
{
  return write_XYZ(stream, point_set);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_POINT_SET_IO_XYZ_H
