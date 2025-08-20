// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef LCC_GEOMETRICAL_TESTS_H
#define LCC_GEOMETRICAL_TESTS_H
///////////////////////////////////////////////////////////////////////////////
#include <CGAL/Linear_cell_complex_operations.h>

namespace lcc_tests
{
///////////////////////////////////////////////////////////////////////////////
enum FACE_ORIENTATION
  {
   ALMOST_XY,
   ALMOST_XZ,
   ALMOST_YZ,
   ORIENTATION_UNKNOWN
  };
const double LCC_GEOMETRICAL_TESTS_EPSILON_DIST=0.000001;
const double LCC_GEOMETRICAL_TESTS_EPSILON_ANGLE=0.1;
///////////////////////////////////////////////////////////////////////////////
bool almost_zero(double d, double epsilon=LCC_GEOMETRICAL_TESTS_EPSILON_DIST)
{ return d>=-epsilon && d<=epsilon; }
///////////////////////////////////////////////////////////////////////////////
bool almost_non_zero(double d, double epsilon=LCC_GEOMETRICAL_TESTS_EPSILON_DIST)
{ return d<-epsilon || d>epsilon; }
///////////////////////////////////////////////////////////////////////////////
bool almost_flat(double angle, double epsilon=LCC_GEOMETRICAL_TESTS_EPSILON_ANGLE)
{ return (angle<epsilon || angle>(360-epsilon) ||
          (angle>180-epsilon && angle<180+epsilon));
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
FACE_ORIENTATION get_face_orientation(LCC& lcc,
                                      typename LCC::Dart_handle dh,
                                      double epsilon=LCC_GEOMETRICAL_TESTS_EPSILON_DIST)
{
  typename LCC::Vector n(CGAL::compute_normal_of_cell_2(lcc, dh));
  if(almost_zero(n.x(), epsilon) && almost_zero(n.y(), epsilon) && almost_non_zero(n.z(), epsilon))
  { return ALMOST_XY; }
  if(almost_zero(n.x(), epsilon) && almost_non_zero(n.y(), epsilon) && almost_zero(n.z(), epsilon))
  { return ALMOST_XZ; }
  if(almost_non_zero(n.x(), epsilon) && almost_zero(n.y(), epsilon) && almost_zero(n.z(), epsilon))
  { return ALMOST_YZ; }
  return ORIENTATION_UNKNOWN;
}
///////////////////////////////////////////////////////////////////////////////
template<typename Point>
bool aligned(const Point& p1, const Point& p2, const Point& p3,
             double epsilon=LCC_GEOMETRICAL_TESTS_EPSILON_ANGLE)
{
  double angle=CGAL::approximate_angle(p1, p2, p3);
  return almost_flat(angle, epsilon);
}
///////////////////////////////////////////////////////////////////////////////
template<typename Vector>
bool coplanar(const Vector& normal1, const Vector& normal2,
             double epsilon=LCC_GEOMETRICAL_TESTS_EPSILON_ANGLE)
{
  double angle=CGAL::approximate_angle(normal1, normal2);
  return almost_flat(angle, epsilon);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
CGAL::Bbox_3 compute_bbox_of_lcc(LCC& pattern)
{
  CGAL::Bbox_3 bbox=pattern.point(pattern.darts().begin()).bbox();
  for(auto it=pattern.darts().begin(), itend=pattern.darts().end(); it!=itend; ++it)
  { bbox+=pattern.point(it).bbox(); }
  return bbox;
}
///////////////////////////////////////////////////////////////////////////////
} // namespace lcc_tests

#endif // LCC_GEOMETRICAL_TESTS_H
