// Copyright (c) 2016 GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_IO_LAS_H
#define CGAL_POINT_SET_IO_LAS_H

#include <CGAL/license/Point_set_3.h>

#ifdef CGAL_LINKED_WITH_LASLIB
#include <CGAL/IO/read_las_points.h>
#include <CGAL/IO/write_las_points.h>
#endif // LAS

#include <fstream>
#include <string>

#if defined(CGAL_LINKED_WITH_LASLIB) || defined(DOXYGEN_RUNNING)

namespace CGAL {

template <typename Point, typename Vector>
class Point_set_3;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {

namespace internal {

template <typename PointSet, typename PropertyMap>
void check_if_property_is_used(PointSet& point_set,
                               PropertyMap& map)
{
  for(typename PointSet::iterator it = point_set.begin(); it != point_set.end(); ++it)
    if(get(map, *it) != typename PropertyMap::value_type())
      return;

  point_set.remove_property_map(map);
}

} // namespace internal

/*!
  \ingroup PkgPointSet3IOLAS

  \brief reads the content of an input stream in the \ref IOStreamLAS into a point set.

  \attention To read a binary file, the flag `std::ios::binary` must be set during the creation of the `ifstream`.

  \param is the input stream
  \param point_set the point set

  \note All LAS properties are read as described in `read_LAS_with_properties()`.

  \return `true` if the reading was successful, `false` otherwise.
 */
template <typename Point, typename Vector>
bool read_LAS(std::istream& is,
              CGAL::Point_set_3<Point, Vector>& point_set)
{
  if(!is)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  typedef CGAL::Point_set_3<Point, Vector> Point_set;
  typedef typename Point_set::template Property_map<float> Float_map;
  typedef typename Point_set::template Property_map<double> Double_map;
  typedef typename Point_set::template Property_map<unsigned short> Ushort_map;
  typedef typename Point_set::template Property_map<unsigned char> Uchar_map;
  typedef typename Point_set::template Property_map<unsigned int> Uint_map;

  Ushort_map intensity = point_set.template add_property_map<unsigned short>("intensity", 0).first;
  Uchar_map return_number = point_set.template add_property_map<unsigned char>("return_number", 0).first;
  Uchar_map number_of_returns = point_set.template add_property_map<unsigned char>("number_of_returns", 0).first;
  Uchar_map scan_direction_flag = point_set.template add_property_map<unsigned char>("scan_direction_flag", 0).first;
  Uchar_map edge_of_flight_line = point_set.template add_property_map<unsigned char>("edge_of_flight_line", 0).first;
  Uchar_map classification = point_set.template add_property_map<unsigned char>("classification", 0).first;
  Uchar_map synthetic_flag = point_set.template add_property_map<unsigned char>("synthetic_flag", 0).first;
  Uchar_map keypoint_flag = point_set.template add_property_map<unsigned char>("keypoint_flag", 0).first;
  Uchar_map withheld_flag = point_set.template add_property_map<unsigned char>("withheld_flag", 0).first;
  Float_map scan_angle = point_set.template add_property_map<float>("scan_angle", 0.).first;
  Uchar_map user_data = point_set.template add_property_map<unsigned char>("user_data", 0).first;
  Ushort_map point_source_ID = point_set.template add_property_map<unsigned short>("point_source_ID", 0).first;
  Uint_map deleted_flag = point_set.template add_property_map<unsigned int>("deleted_flag", 0).first;
  Double_map gps_time = point_set.template add_property_map<double>("gps_time", 0).first;
  Ushort_map R = point_set.template add_property_map<unsigned short>("R", 0).first;
  Ushort_map G = point_set.template add_property_map<unsigned short>("G", 0).first;
  Ushort_map B = point_set.template add_property_map<unsigned short>("B", 0).first;
  Ushort_map I = point_set.template add_property_map<unsigned short>("I", 0).first;

  bool okay
      = read_LAS_with_properties
      (is, point_set.index_back_inserter(),
       make_las_point_reader(point_set.point_push_map()),
       std::make_pair(point_set.push_property_map(intensity), LAS_property::Intensity()),
       std::make_pair(point_set.push_property_map(return_number), LAS_property::Return_number()),
       std::make_pair(point_set.push_property_map(number_of_returns), LAS_property::Number_of_returns()),
       std::make_pair(point_set.push_property_map(scan_direction_flag), LAS_property::Scan_direction_flag()),
       std::make_pair(point_set.push_property_map(edge_of_flight_line), LAS_property::Edge_of_flight_line()),
       std::make_pair(point_set.push_property_map(classification), LAS_property::Classification()),
       std::make_pair(point_set.push_property_map(synthetic_flag), LAS_property::Synthetic_flag()),
       std::make_pair(point_set.push_property_map(keypoint_flag), LAS_property::Keypoint_flag()),
       std::make_pair(point_set.push_property_map(withheld_flag), LAS_property::Withheld_flag()),
       std::make_pair(point_set.push_property_map(scan_angle), LAS_property::Scan_angle()),
       std::make_pair(point_set.push_property_map(user_data), LAS_property::User_data()),
       std::make_pair(point_set.push_property_map(point_source_ID), LAS_property::Point_source_ID()),
       std::make_pair(point_set.push_property_map(deleted_flag), LAS_property::Deleted_flag()),
       std::make_pair(point_set.push_property_map(gps_time), LAS_property::GPS_time()),
       std::make_pair(point_set.push_property_map(R), LAS_property::R()),
       std::make_pair(point_set.push_property_map(G), LAS_property::G()),
       std::make_pair(point_set.push_property_map(B), LAS_property::B()),
       std::make_pair(point_set.push_property_map(I), LAS_property::I()));

  internal::check_if_property_is_used(point_set, intensity);
  internal::check_if_property_is_used(point_set, return_number);
  internal::check_if_property_is_used(point_set, number_of_returns);
  internal::check_if_property_is_used(point_set, scan_direction_flag);
  internal::check_if_property_is_used(point_set, edge_of_flight_line);
  internal::check_if_property_is_used(point_set, classification);
  internal::check_if_property_is_used(point_set, synthetic_flag);
  internal::check_if_property_is_used(point_set, keypoint_flag);
  internal::check_if_property_is_used(point_set, withheld_flag);
  internal::check_if_property_is_used(point_set, scan_angle);
  internal::check_if_property_is_used(point_set, user_data);
  internal::check_if_property_is_used(point_set, point_source_ID);
  internal::check_if_property_is_used(point_set, deleted_flag);
  internal::check_if_property_is_used(point_set, gps_time);
  internal::check_if_property_is_used(point_set, R);
  internal::check_if_property_is_used(point_set, G);
  internal::check_if_property_is_used(point_set, B);
  internal::check_if_property_is_used(point_set, I);

  return okay;
}

/*!
  \ingroup PkgPointSet3IOLAS

  \brief reads the content of an input file in the \ref IOStreamLAS into a point set.

  \param fname the path to the input file
  \param point_set the point set

  \note All LAS properties are read as described in `read_LAS_with_properties()`.

  \return `true` if the reading was successful, `false` otherwise.
*/
template <typename Point, typename Vector>
bool read_LAS(const std::string& fname, CGAL::Point_set_3<Point, Vector>& point_set)
{
  std::ifstream is(fname, std::ios::binary);
  CGAL::IO::set_mode(is, CGAL::IO::BINARY);
  return read_LAS(is, point_set);
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

template <typename Point, typename Vector>
CGAL_DEPRECATED bool read_las_point_set(std::istream& is, ///< input stream.
                                        CGAL::Point_set_3<Point, Vector>& point_set) ///< point set
{
  return IO::read_LAS(is, point_set);
}

#endif // CGAL_NO_DEPRECATED_CODE

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace IO {

/*!
  \ingroup PkgPointSet3IOLAS

  \brief writes the content of a point set into an output stream in the \ref IOStreamLAS.

  \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation of the `ofstream`.

  \tparam Point the point type of the `Point_set_3`
  \tparam Vector the vector type of the `Point_set_3`

  \param os the output stream
  \param point_set the point set

  \note All LAS properties are written as described in `read_LAS_with_properties()`.

  \return `true` if the writing was successful, `false` otherwise.
 */
template <typename Point, typename Vector>
bool write_LAS(std::ostream& os,
               CGAL::Point_set_3<Point, Vector>& point_set)
{
  if(!os)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  typedef CGAL::Point_set_3<Point, Vector> Point_set;
  typedef typename Point_set::template Property_map<float> Float_map;
  typedef typename Point_set::template Property_map<double> Double_map;
  typedef typename Point_set::template Property_map<unsigned short> Ushort_map;
  typedef typename Point_set::template Property_map<unsigned char> Uchar_map;
  typedef typename Point_set::template Property_map<unsigned int> Uint_map;

  Ushort_map intensity;
  bool remove_intensity;
  boost::tie(intensity, remove_intensity)
      = point_set.template add_property_map<unsigned short>("intensity", 0);

  Uchar_map return_number;
  bool remove_return_number;
  boost::tie(return_number, remove_return_number)
      = point_set.template add_property_map<unsigned char>("return_number", 0);

  Uchar_map number_of_returns;
  bool remove_number_of_returns;
  boost::tie(number_of_returns, remove_number_of_returns)
      = point_set.template add_property_map<unsigned char>("number_of_returns", 0);

  Uchar_map scan_direction_flag;
  bool remove_scan_direction_flag;
  boost::tie(scan_direction_flag, remove_scan_direction_flag)
      = point_set.template add_property_map<unsigned char>("scan_direction_flag", 0);

  Uchar_map edge_of_flight_line;
  bool remove_edge_of_flight_line;
  boost::tie(edge_of_flight_line, remove_edge_of_flight_line)
      = point_set.template add_property_map<unsigned char>("edge_of_flight_line", 0);

  Uchar_map classification;
  bool remove_classification;
  boost::tie(classification, remove_classification)
      = point_set.template add_property_map<unsigned char>("classification", 0);

  Uchar_map synthetic_flag;
  bool remove_synthetic_flag;
  boost::tie(synthetic_flag, remove_synthetic_flag)
      = point_set.template add_property_map<unsigned char>("synthetic_flag", 0);

  Uchar_map keypoint_flag;
  bool remove_keypoint_flag;
  boost::tie(keypoint_flag, remove_keypoint_flag)
      = point_set.template add_property_map<unsigned char>("keypoint_flag", 0);

  Uchar_map withheld_flag;
  bool remove_withheld_flag;
  boost::tie(withheld_flag, remove_withheld_flag)
      = point_set.template add_property_map<unsigned char>("withheld_flag", 0);

  Float_map scan_angle;
  bool remove_scan_angle;
  boost::tie(scan_angle, remove_scan_angle)
      = point_set.template add_property_map<float>("scan_angle", 0.);

  Uchar_map user_data;
  bool remove_user_data;
  boost::tie(user_data, remove_user_data)
      = point_set.template add_property_map<unsigned char>("user_data", 0);

  Ushort_map point_source_ID;
  bool remove_point_source_ID;
  boost::tie(point_source_ID, remove_point_source_ID)
      = point_set.template add_property_map<unsigned short>("point_source_ID", 0);

  Uint_map deleted_flag;
  bool remove_deleted_flag;
  boost::tie(deleted_flag, remove_deleted_flag)
      = point_set.template add_property_map<unsigned int>("deleted_flag", 0);

  Double_map gps_time;
  bool remove_gps_time;
  boost::tie(gps_time, remove_gps_time)
      = point_set.template add_property_map<double>("gps_time", 0);

  Ushort_map R;
  bool remove_R;
  boost::tie(R, remove_R) = point_set.template add_property_map<unsigned short>("R", 0);
  Ushort_map G;
  bool remove_G;
  boost::tie(G, remove_G) = point_set.template add_property_map<unsigned short>("G", 0);
  Ushort_map B;
  bool remove_B;
  boost::tie(B, remove_B) = point_set.template add_property_map<unsigned short>("B", 0);
  Ushort_map I;
  bool remove_I;
  boost::tie(I, remove_I) = point_set.template add_property_map<unsigned short>("I", 0);

  if(remove_R)
  {
    Uchar_map charR, charG, charB;
    bool foundR, foundG, foundB;
    boost::tie(charR, foundR) = point_set.template property_map<unsigned char>("r");
    if(!foundR)
      boost::tie(charR, foundR) = point_set.template property_map<unsigned char>("red");
    boost::tie(charG, foundG) = point_set.template property_map<unsigned char>("g");
    if(!foundG)
      boost::tie(charG, foundG) = point_set.template property_map<unsigned char>("green");
    boost::tie(charB, foundB) = point_set.template property_map<unsigned char>("b");
    if(!foundB)
      boost::tie(charB, foundB) = point_set.template property_map<unsigned char>("blue");

    if(foundR && foundG && foundB)
    {
      for(typename Point_set::iterator it = point_set.begin(); it != point_set.end(); ++it)
      {
        put(R, *it, (unsigned short)(get(charR, *it)));
        put(G, *it, (unsigned short)(get(charG, *it)));
        put(B, *it, (unsigned short)(get(charB, *it)));
      }
    }
  }

  bool okay
      = write_LAS_with_properties
      (os, point_set,
       make_las_point_writer(point_set.point_map()),
       std::make_pair(intensity, LAS_property::Intensity()),
       std::make_pair(return_number, LAS_property::Return_number()),
       std::make_pair(number_of_returns, LAS_property::Number_of_returns()),
       std::make_pair(scan_direction_flag, LAS_property::Scan_direction_flag()),
       std::make_pair(edge_of_flight_line, LAS_property::Edge_of_flight_line()),
       std::make_pair(classification, LAS_property::Classification()),
       std::make_pair(synthetic_flag, LAS_property::Synthetic_flag()),
       std::make_pair(keypoint_flag, LAS_property::Keypoint_flag()),
       std::make_pair(withheld_flag, LAS_property::Withheld_flag()),
       std::make_pair(scan_angle, LAS_property::Scan_angle()),
       std::make_pair(user_data, LAS_property::User_data()),
       std::make_pair(point_source_ID, LAS_property::Point_source_ID()),
       std::make_pair(deleted_flag, LAS_property::Deleted_flag()),
       std::make_pair(gps_time, LAS_property::GPS_time()),
       std::make_pair(R, LAS_property::R()),
       std::make_pair(G, LAS_property::G()),
       std::make_pair(B, LAS_property::B()),
       std::make_pair(I, LAS_property::I()));

  if(remove_intensity) point_set.remove_property_map(intensity);
  if(remove_return_number) point_set.remove_property_map(return_number);
  if(remove_number_of_returns) point_set.remove_property_map(number_of_returns);
  if(remove_scan_direction_flag) point_set.remove_property_map(scan_direction_flag);
  if(remove_edge_of_flight_line) point_set.remove_property_map(edge_of_flight_line);
  if(remove_classification) point_set.remove_property_map(classification);
  if(remove_synthetic_flag) point_set.remove_property_map(synthetic_flag);
  if(remove_keypoint_flag) point_set.remove_property_map(keypoint_flag);
  if(remove_withheld_flag) point_set.remove_property_map(withheld_flag);
  if(remove_scan_angle) point_set.remove_property_map(scan_angle);
  if(remove_user_data) point_set.remove_property_map(user_data);
  if(remove_point_source_ID) point_set.remove_property_map(point_source_ID);
  if(remove_deleted_flag) point_set.remove_property_map(deleted_flag);
  if(remove_gps_time) point_set.remove_property_map(gps_time);
  if(remove_R) point_set.remove_property_map(R);
  if(remove_G) point_set.remove_property_map(G);
  if(remove_B) point_set.remove_property_map(B);
  if(remove_I) point_set.remove_property_map(I);

  return okay;
}

/*!
  \ingroup PkgPointSet3IOLAS

  \brief writes the content of a point set into an output file in the \ref IOStreamLAS.

  \tparam Point the point type of the `Point_set_3`
  \tparam Vector the vector type of the `Point_set_3`

  \param fname the path to the output file
  \param point_set the point set

  \note All LAS properties are written as described in `read_LAS_with_properties()`.

  \return `true` if the writing was successful, `false` otherwise.
 */
template <typename Point, typename Vector>
bool write_LAS(const std::string& fname,
               CGAL::Point_set_3<Point, Vector>& point_set)
{
  std::ofstream os(fname, std::ios::binary);
  CGAL::IO::set_mode(os, CGAL::IO::BINARY);
  return write_LAS(os, point_set);
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

template <typename Point, typename Vector>
CGAL_DEPRECATED bool write_las_point_set(std::ostream& os, ///< output stream.
                                         CGAL::Point_set_3<Point, Vector>& point_set)  ///< point set
{
  return IO::write_LAS(os, point_set);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // LAS

#endif // CGAL_POINT_SET_IO_LAS_H
