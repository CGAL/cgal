// Copyright (c) 2016  GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_3_IO
#define CGAL_POINT_SET_3_IO

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/read_ply_point_set_3.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/write_ply_points.h>


namespace CGAL {

/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point, typename Vector>
bool
read_xyz_point_set(
  std::istream& stream, ///< input stream.
  Point_set_3<Point, Vector>& point_set) ///< point set
{
  point_set.add_normal_map();

  bool out = CGAL::read_xyz_points_and_normals
    (stream,
     point_set.index_back_inserter(),
     point_set.point_push_map(),
     point_set.normal_push_map());

  bool has_normals = false;
  for (typename Point_set_3<Point, Vector>::const_iterator it = point_set.begin();
       it != point_set.end(); ++ it)
    if (point_set.normal(*it) != CGAL::NULL_VECTOR)
      {
        has_normals = true;
        break;
      }

  if (!has_normals)
    point_set.remove_normal_map();
  
  return out;
}

/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point, typename Vector>
bool
read_off_point_set(
  std::istream& stream, ///< input stream.
  Point_set_3<Point, Vector>& point_set) ///< point set
{
  point_set.add_normal_map();

  bool out = CGAL::read_off_points_and_normals
    (stream,
     point_set.index_back_inserter(),
     point_set.point_push_map(),
     point_set.normal_push_map());

  bool has_normals = false;
  for (typename Point_set_3<Point, Vector>::const_iterator it = point_set.begin();
       it != point_set.end(); ++ it)
    if (point_set.normal(*it) != CGAL::NULL_VECTOR)
      {
        has_normals = true;
        break;
      }

  if (!has_normals)
    point_set.remove_normal_map();

  return out;
}


  
/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point, typename Vector>
bool
read_ply_point_set(
  std::istream& stream, ///< input stream.
  Point_set_3<Point, Vector>& point_set) ///< point set
{
  CGAL::Ply_interpreter_point_set_3<Point, Vector> interpreter (point_set);

  return CGAL::read_ply_custom_points
    (stream, interpreter,
     typename Kernel_traits<Point>::Kernel());
}

/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point, typename Vector>
bool
write_xyz_point_set(
  std::ostream& stream, ///< output stream.
  const Point_set_3<Point, Vector>& point_set)  ///< point set
{
  if (point_set.has_normal_map())
    return CGAL::write_xyz_points_and_normals
      (stream, point_set.begin(), point_set.end(),
       point_set.point_map(), point_set.normal_map());
  
  return CGAL::write_xyz_points
  (stream, point_set.begin(), point_set.end(),
   point_set.point_map());
}

/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point, typename Vector>
bool
write_off_point_set(
  std::ostream& stream, ///< output stream.
  const Point_set_3<Point, Vector>& point_set)  ///< point set
{
  if (point_set.has_normal_map())
    return CGAL::write_off_points_and_normals
      (stream, point_set.begin(), point_set.end(),
       point_set.point_map(), point_set.normal_map());
  
  return CGAL::write_off_points
  (stream, point_set.begin(), point_set.end(),
   point_set.point_map());
}

/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point, typename Vector>
bool
write_ply_point_set(
  std::ostream& stream, ///< output stream.
  const Point_set_3<Point, Vector>& point_set)  ///< point set
{

  stream << point_set;
  return true;
  // if (point_set.has_normal_map())
  //   return CGAL::write_ply_points_and_normals
  //     (stream, point_set.begin(), point_set.end(),
  //      point_set.point_map(), point_set.normal_map());
  
  // return CGAL::write_ply_points
  // (stream, point_set.begin(), point_set.end(),
  //  point_set.point_map());
}

} // namespace CGAL


#endif // CGAL_POINT_SET_3_IO
