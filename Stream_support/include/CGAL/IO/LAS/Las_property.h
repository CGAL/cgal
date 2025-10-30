// Copyright (c) 2017  GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_IO_LAS_LAS_PROPERTY_H
#define CGAL_IO_LAS_LAS_PROPERTY_H

namespace CGAL {
  namespace IO {

namespace LAS_property {
namespace Id {

enum Id
{
  X,
  Y,
  Z,
  Intensity,
  Return_number,
  Number_of_returns,
  Scan_direction_flag,
  Edge_of_flight_line,
  Classification,
  Synthetic_flag,
  Keypoint_flag,
  Withheld_flag,
  Scan_angle,
  User_data,
  Point_source_ID,
  Deleted_flag,
  GPS_time,
  R,
  G,
  B,
  I
};

} // namespace Id

template <typename T, Id::Id id>
struct Base
{
  typedef T type;
};

typedef Base<double, Id::X> X;
typedef Base<double, Id::Y> Y;
typedef Base<double, Id::Z> Z;
typedef Base<unsigned short, Id::Intensity> Intensity;
typedef Base<unsigned char, Id::Return_number> Return_number;
typedef Base<unsigned char, Id::Number_of_returns> Number_of_returns;
typedef Base<unsigned char, Id::Scan_direction_flag> Scan_direction_flag;
typedef Base<unsigned char, Id::Edge_of_flight_line> Edge_of_flight_line;
typedef Base<unsigned char, Id::Classification> Classification;
typedef Base<unsigned char, Id::Synthetic_flag> Synthetic_flag;
typedef Base<unsigned char, Id::Keypoint_flag> Keypoint_flag;
typedef Base<unsigned char, Id::Withheld_flag> Withheld_flag;
typedef Base<float, Id::Scan_angle> Scan_angle;
typedef Base<unsigned char, Id::User_data> User_data;
typedef Base<unsigned short, Id::Point_source_ID> Point_source_ID;
typedef Base<unsigned int, Id::Deleted_flag> Deleted_flag;
typedef Base<double, Id::GPS_time> GPS_time;
typedef Base<unsigned short, Id::R> R;
typedef Base<unsigned short, Id::G> G;
typedef Base<unsigned short, Id::B> B;
typedef Base<unsigned short, Id::I> I;
} // namespace LAS_property

// documenation in ../LAS.h
template <typename PointMap>
std::tuple<PointMap,
           typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3,
           LAS_property::X, LAS_property::Y, LAS_property::Z >
make_las_point_reader(PointMap point_map)
{
  return std::make_tuple (point_map, typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3(),
                          LAS_property::X(), LAS_property::Y(), LAS_property::Z());
}


} // namespace IO
} // namespace CGAL

#endif
