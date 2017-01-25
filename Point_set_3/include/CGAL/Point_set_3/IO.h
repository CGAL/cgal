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

#include <CGAL/license/Point_set_3.h>


#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/read_ply_point_set_3.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/write_ply_points.h>


#ifdef CGAL_LINKED_WITH_LASLIB
#include <CGAL/IO/read_las_points.h>
#include <CGAL/IO/write_las_points.h>
#endif

namespace CGAL {

/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point, typename Vector>
bool
read_xyz_point_set(
  std::istream& stream, ///< input stream.
  CGAL::Point_set_3<Point, Vector>& point_set) ///< point set
{
  point_set.add_normal_map();

  bool out = CGAL::read_xyz_points_and_normals
    (stream,
     point_set.index_back_inserter(),
     point_set.point_push_map(),
     point_set.normal_push_map());

  bool has_normals = false;
  for (typename CGAL::Point_set_3<Point, Vector>::const_iterator it = point_set.begin();
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
  CGAL::Point_set_3<Point, Vector>& point_set) ///< point set
{
  point_set.add_normal_map();

  bool out = CGAL::read_off_points_and_normals
    (stream,
     point_set.index_back_inserter(),
     point_set.point_push_map(),
     point_set.normal_push_map());

  bool has_normals = false;
  for (typename CGAL::Point_set_3<Point, Vector>::const_iterator it = point_set.begin();
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
  CGAL::Point_set_3<Point, Vector>& point_set) ///< point set
{
  if(!stream)
    {
      std::cerr << "Error: cannot open file" << std::endl;
      return false;
    }

  Ply::internal::Ply_reader reader;
  
  if (!(reader.init (stream)))
    return false;
  
  Ply::internal::Ply_interpreter_point_set_3<Point, Vector> interpreter (point_set);
  interpreter.instantiate_properties (reader);
  
  std::size_t points_read = 0;
  
  while (!(stream.eof()) && points_read < reader.m_nb_points)
    {
      for (std::size_t i = 0; i < reader.readers().size (); ++ i)
        reader.readers()[i]->get (stream);

      interpreter.process_line (reader);
      
      ++ points_read;
    }

  return (points_read == reader.m_nb_points);
}

/*!
  \ingroup PkgPointSet3IO
 */
#ifdef CGAL_LINKED_WITH_LASLIB

namespace internal
{
  template <typename PointSet, typename PropertyMap>
  void check_if_property_is_used (PointSet& point_set,
                                  PropertyMap& map)
  {
    for (typename PointSet::iterator it = point_set.begin(); it != point_set.end(); ++ it)
      if (get(map, *it) != typename PropertyMap::value_type())
        return;

    point_set.remove_property_map (map);
  }

}
  
template <typename Point, typename Vector>
bool
read_las_point_set(
  std::istream& stream, ///< input stream.
  CGAL::Point_set_3<Point, Vector>& point_set) ///< point set
{
  if(!stream)
    {
      std::cerr << "Error: cannot open file" << std::endl;
      return false;
    }

  typedef CGAL::Point_set_3<Point, Vector> Point_set;
  typedef typename Point_set::template Property_map<float> Float_map;
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
  Ushort_map R = point_set.template add_property_map<unsigned short>("R", 0).first;
  Ushort_map G = point_set.template add_property_map<unsigned short>("G", 0).first;
  Ushort_map B = point_set.template add_property_map<unsigned short>("B", 0).first;
  Ushort_map I = point_set.template add_property_map<unsigned short>("I", 0).first;

  bool okay
    = read_las_points_with_properties
    (stream, point_set.index_back_inserter(),
     Las::point_property (point_set.point_push_map()),
     std::make_pair (point_set.push_property_map (intensity), Las::Property::intensity()),
     std::make_pair (point_set.push_property_map (return_number), Las::Property::return_number()),
     std::make_pair (point_set.push_property_map (number_of_returns), Las::Property::number_of_returns()),
     std::make_pair (point_set.push_property_map (scan_direction_flag), Las::Property::scan_direction_flag()),
     std::make_pair (point_set.push_property_map (edge_of_flight_line), Las::Property::edge_of_flight_line()),
     std::make_pair (point_set.push_property_map (classification), Las::Property::classification()),
     std::make_pair (point_set.push_property_map (synthetic_flag), Las::Property::synthetic_flag()),
     std::make_pair (point_set.push_property_map (keypoint_flag), Las::Property::keypoint_flag()),
     std::make_pair (point_set.push_property_map (withheld_flag), Las::Property::withheld_flag()),
     std::make_pair (point_set.push_property_map (scan_angle), Las::Property::scan_angle()),
     std::make_pair (point_set.push_property_map (user_data), Las::Property::user_data()),
     std::make_pair (point_set.push_property_map (point_source_ID), Las::Property::point_source_ID()),
     std::make_pair (point_set.push_property_map (deleted_flag), Las::Property::deleted_flag()),
     std::make_pair (point_set.push_property_map (R), Las::Property::R()),
     std::make_pair (point_set.push_property_map (G), Las::Property::G()),
     std::make_pair (point_set.push_property_map (B), Las::Property::B()),
     std::make_pair (point_set.push_property_map (I), Las::Property::I()));

  internal::check_if_property_is_used (point_set, intensity);
  internal::check_if_property_is_used (point_set, return_number);
  internal::check_if_property_is_used (point_set, number_of_returns);
  internal::check_if_property_is_used (point_set, scan_direction_flag);
  internal::check_if_property_is_used (point_set, edge_of_flight_line);
  internal::check_if_property_is_used (point_set, classification);
  internal::check_if_property_is_used (point_set, synthetic_flag);
  internal::check_if_property_is_used (point_set, keypoint_flag);
  internal::check_if_property_is_used (point_set, withheld_flag);
  internal::check_if_property_is_used (point_set, scan_angle);
  internal::check_if_property_is_used (point_set, user_data);
  internal::check_if_property_is_used (point_set, point_source_ID);
  internal::check_if_property_is_used (point_set, deleted_flag);
  internal::check_if_property_is_used (point_set, R);
  internal::check_if_property_is_used (point_set, G);
  internal::check_if_property_is_used (point_set, B);
  internal::check_if_property_is_used (point_set, I);
  
  return okay;
}

/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point, typename Vector>
bool
write_las_point_set(
  std::ostream& stream, ///< output stream.
  CGAL::Point_set_3<Point, Vector>& point_set)  ///< point set
{
  if(!stream)
    {
      std::cerr << "Error: cannot open file" << std::endl;
      return false;
    }

  typedef CGAL::Point_set_3<Point, Vector> Point_set;
  typedef typename Point_set::template Property_map<float> Float_map;
  typedef typename Point_set::template Property_map<unsigned short> Ushort_map;
  typedef typename Point_set::template Property_map<unsigned char> Uchar_map;
  typedef typename Point_set::template Property_map<unsigned int> Uint_map;

  Ushort_map intensity;
  bool remove_intensity;
  boost::tie(intensity, remove_intensity)
    = point_set.template add_property_map<unsigned short>("intensity", 0);
  
  Uchar_map return_number;
  bool remove_return_number;
  boost::tie (return_number, remove_return_number)
    = point_set.template add_property_map<unsigned char>("return_number", 0);
  
  Uchar_map number_of_returns;
  bool remove_number_of_returns;
  boost::tie (number_of_returns, remove_number_of_returns)
    = point_set.template add_property_map<unsigned char>("number_of_returns", 0);
    
  Uchar_map scan_direction_flag;
  bool remove_scan_direction_flag;
  boost::tie (scan_direction_flag, remove_scan_direction_flag)
    = point_set.template add_property_map<unsigned char>("scan_direction_flag", 0);
  
  Uchar_map edge_of_flight_line;
  bool remove_edge_of_flight_line;
  boost::tie (edge_of_flight_line, remove_edge_of_flight_line)
    = point_set.template add_property_map<unsigned char>("edge_of_flight_line", 0);
  
  Uchar_map classification;
  bool remove_classification;
  boost::tie (classification, remove_classification)
    = point_set.template add_property_map<unsigned char>("classification", 0);
  
  Uchar_map synthetic_flag;
  bool remove_synthetic_flag;
  boost::tie (synthetic_flag, remove_synthetic_flag)
    = point_set.template add_property_map<unsigned char>("synthetic_flag", 0);
  
  Uchar_map keypoint_flag;
  bool remove_keypoint_flag;
  boost::tie (keypoint_flag, remove_keypoint_flag)
    = point_set.template add_property_map<unsigned char>("keypoint_flag", 0);
  
  Uchar_map withheld_flag;
  bool remove_withheld_flag;
  boost::tie (withheld_flag, remove_withheld_flag)
    = point_set.template add_property_map<unsigned char>("withheld_flag", 0);
  
  Float_map scan_angle;
  bool remove_scan_angle;
  boost::tie (scan_angle, remove_scan_angle)
    = point_set.template add_property_map<float>("scan_angle", 0.);
  
  Uchar_map user_data;
  bool remove_user_data;
  boost::tie (user_data, remove_user_data)
    = point_set.template add_property_map<unsigned char>("user_data", 0);
  
  Ushort_map point_source_ID;
  bool remove_point_source_ID;
  boost::tie (point_source_ID, remove_point_source_ID)
    = point_set.template add_property_map<unsigned short>("point_source_ID", 0);
  
  Uint_map deleted_flag;
  bool remove_deleted_flag;
  boost::tie (deleted_flag, remove_deleted_flag)
    = point_set.template add_property_map<unsigned int>("deleted_flag", 0);
  
  Ushort_map R;
  bool remove_R;
  boost::tie (R, remove_R) = point_set.template add_property_map<unsigned short>("R", 0);
  Ushort_map G;
  bool remove_G;
  boost::tie (G, remove_G) = point_set.template add_property_map<unsigned short>("G", 0);
  Ushort_map B;
  bool remove_B;
  boost::tie (B, remove_B) = point_set.template add_property_map<unsigned short>("B", 0);
  Ushort_map I;
  bool remove_I;
  boost::tie (I, remove_I) = point_set.template add_property_map<unsigned short>("I", 0);

  if (remove_R)
    {
      Uchar_map charR, charG, charB;
      bool foundR, foundG, foundB;
      boost::tie (charR, foundR) = point_set.template property_map<unsigned char>("r");
      if (!foundR)
        boost::tie (charR, foundR) = point_set.template property_map<unsigned char>("red");
      boost::tie (charG, foundG) = point_set.template property_map<unsigned char>("g");
      if (!foundG)
        boost::tie (charG, foundG) = point_set.template property_map<unsigned char>("green");
      boost::tie (charB, foundB) = point_set.template property_map<unsigned char>("b");
      if (!foundB)
        boost::tie (charB, foundB) = point_set.template property_map<unsigned char>("blue");

      if (foundR && foundG && foundB)
        {
          for (typename Point_set::iterator it = point_set.begin(); it != point_set.end(); ++ it)
            {
              put (R, *it, (unsigned short)(get(charR, *it)));
              put (G, *it, (unsigned short)(get(charG, *it)));
              put (B, *it, (unsigned short)(get(charB, *it)));
            }
        }
    }
  
  bool okay
    = write_las_points_with_properties
    (stream, point_set.begin(), point_set.end(),
     Las::point_writer (point_set.point_map()),
     std::make_pair (intensity, Las::Property::intensity()),
     std::make_pair (return_number, Las::Property::return_number()),
     std::make_pair (number_of_returns, Las::Property::number_of_returns()),
     std::make_pair (scan_direction_flag, Las::Property::scan_direction_flag()),
     std::make_pair (edge_of_flight_line, Las::Property::edge_of_flight_line()),
     std::make_pair (classification, Las::Property::classification()),
     std::make_pair (synthetic_flag, Las::Property::synthetic_flag()),
     std::make_pair (keypoint_flag, Las::Property::keypoint_flag()),
     std::make_pair (withheld_flag, Las::Property::withheld_flag()),
     std::make_pair (scan_angle, Las::Property::scan_angle()),
     std::make_pair (user_data, Las::Property::user_data()),
     std::make_pair (point_source_ID, Las::Property::point_source_ID()),
     std::make_pair (deleted_flag, Las::Property::deleted_flag()),
     std::make_pair (R, Las::Property::R()),
     std::make_pair (G, Las::Property::G()),
     std::make_pair (B, Las::Property::B()),
     std::make_pair (I, Las::Property::I()));

  if (remove_intensity) point_set.remove_property_map (intensity);
  if (remove_return_number) point_set.remove_property_map (return_number);
  if (remove_number_of_returns) point_set.remove_property_map (number_of_returns);
  if (remove_scan_direction_flag) point_set.remove_property_map (scan_direction_flag);
  if (remove_edge_of_flight_line) point_set.remove_property_map (edge_of_flight_line);
  if (remove_classification) point_set.remove_property_map (classification);
  if (remove_synthetic_flag) point_set.remove_property_map (synthetic_flag);
  if (remove_keypoint_flag) point_set.remove_property_map (keypoint_flag);
  if (remove_withheld_flag) point_set.remove_property_map (withheld_flag);
  if (remove_scan_angle) point_set.remove_property_map (scan_angle);
  if (remove_user_data) point_set.remove_property_map (user_data);
  if (remove_point_source_ID) point_set.remove_property_map (point_source_ID);
  if (remove_deleted_flag) point_set.remove_property_map (deleted_flag);
  if (remove_R) point_set.remove_property_map (R);
  if (remove_G) point_set.remove_property_map (G);
  if (remove_B) point_set.remove_property_map (B);
  if (remove_I) point_set.remove_property_map (I);
  
  return okay;
}
  
#endif
  
/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point, typename Vector>
bool
write_xyz_point_set(
  std::ostream& stream, ///< output stream.
  const CGAL::Point_set_3<Point, Vector>& point_set)  ///< point set
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
  const CGAL::Point_set_3<Point, Vector>& point_set)  ///< point set
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
  const CGAL::Point_set_3<Point, Vector>& point_set)  ///< point set
{

  stream << point_set;
  return true;
}


/*!
  
  \ingroup PkgPointSet3IO

  \brief Reads the point set from an input stream that can be either:

  - XYZ
  - OFF
  - PLY

  The format is detected from the stream. If the stream contains
  normal vectors, the normal map is added to the point set. For PLY
  input, all point properties found in the header are added.
  \relates Point_set_3
*/
template <typename Point, typename Vector>
std::istream& operator>>(std::istream& is,
                         CGAL::Point_set_3<Point, Vector>& ps)
{
  // Check format identifier on first line
  std::string line;
  if (!getline(is, line))
    return is;

  is.seekg(0);
  if (line.find("OFF") == 0 || line.find("NOFF") == 0)
    CGAL::read_off_point_set (is, ps);
  else if (line.find("ply") == 0)
    CGAL::read_ply_point_set (is, ps);
#ifdef CGAL_LINKED_WITH_LASLIB
  else if (line == "LASF")
    CGAL::read_las_point_set (is, ps);
#endif
  else
    CGAL::read_xyz_point_set (is, ps);
    
  return is;
}

/// \cond SKIP_IN_MANUAL
namespace internal
{
  template <typename Point, typename Vector>
  class Abstract_property_printer
  {
  public:
    virtual ~Abstract_property_printer() { }
    virtual void print (std::ostream& stream, const typename CGAL::Point_set_3<Point,Vector>::Index& index) = 0;
  };

  template <typename Point, typename Vector, typename Type>
  class Property_printer : public Abstract_property_printer<Point, Vector>
  {
    typedef typename CGAL::Point_set_3<Point, Vector> Point_set;
    typedef typename Point_set::template Property_map<Type> Pmap;
    Pmap m_pmap;
  public:
    Property_printer (const Pmap& pmap) : m_pmap (pmap)
    {

    }
    
    virtual void print(std::ostream& stream, const typename CGAL::Point_set_3<Point,Vector>::Index& index)
    {
      stream << get(m_pmap, index);
    }
  };

  template <typename Point, typename Vector, typename Type>
  class Simple_property_printer : public Abstract_property_printer<Point, Vector>
  {
    typedef typename CGAL::Point_set_3<Point, Vector> Point_set;
    typedef typename Point_set::template Property_map<Type> Pmap;
    Pmap m_pmap;
  public:
    Simple_property_printer (const Pmap& pmap) : m_pmap (pmap)
    {

    }
    
    virtual void print(std::ostream& stream, const typename CGAL::Point_set_3<Point,Vector>::Index& index)
    {
      if (get_mode(stream) == IO::ASCII)
        stream << get(m_pmap, index);
      else
        {
          Type t = get (m_pmap, index);
          stream.write (reinterpret_cast<char*>(&t), sizeof(t));
        }
    }
  };

  template <typename Point, typename Vector, typename Type>
  class Char_property_printer : public Abstract_property_printer<Point, Vector>
  {
    typedef typename CGAL::Point_set_3<Point, Vector> Point_set;
    typedef typename Point_set::template Property_map<Type> Pmap;
    Pmap m_pmap;
  public:
    Char_property_printer (const Pmap& pmap) : m_pmap (pmap)
    {

    }
    
    virtual void print(std::ostream& stream, const typename CGAL::Point_set_3<Point,Vector>::Index& index)
    {
      if (get_mode(stream) == IO::ASCII)
        stream << int(get(m_pmap, index));
      else
        {
          Type t = get (m_pmap, index);
          stream.write (reinterpret_cast<char*>(&t), sizeof(t));
        }
    }
  };
  
}
/// \endcond
  
/*!
  
  \ingroup PkgPointSet3IO

  \brief Inserts the point set in an output stream in ASCII PLY
  format. All properties are inserted in their instantiation order.

  \relates Point_set_3
*/
template <typename Point, typename Vector>
std::ostream& operator<<(std::ostream& os,
                         const CGAL::Point_set_3<Point, Vector>& ps)
{
  typedef CGAL::Point_set_3<Point, Vector> Point_set;
    
  os << "ply" << std::endl
     << ((get_mode(os) == IO::BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
     << "comment Generated by the CGAL library" << std::endl
     << "element vertex " << ps.number_of_points() << std::endl;
  
  std::vector<std::string> prop = ps.base().properties();
  std::vector<internal::Abstract_property_printer<Point, Vector>*> printers;
  
  for (std::size_t i = 0; i < prop.size(); ++ i)
    {
      if (prop[i] == "index")
        continue;

      if (prop[i] == "point")
        {
          os << "property double x" << std::endl
             << "property double y" << std::endl
             << "property double z" << std::endl;
          printers.push_back (new internal::Property_printer<Point,Vector,Point>(ps.point_map()));
          continue;
        }
      if (prop[i] == "normal")
        {
          os << "property double nx" << std::endl
             << "property double ny" << std::endl
             << "property double nz" << std::endl;
          printers.push_back (new internal::Property_printer<Point,Vector,Vector>(ps.normal_map()));
          continue;
        }
      
      bool okay = false;
      {
        typename Point_set::template Property_map<boost::int8_t> pmap;
        boost::tie (pmap, okay) = ps.template property_map<boost::int8_t>(prop[i]);
        if (okay)
          {
            os << "property char " << prop[i] << std::endl;
            printers.push_back (new internal::Char_property_printer<Point,Vector,boost::int8_t>(pmap));
            continue;
          }
      }
      {
        typename Point_set::template Property_map<boost::uint8_t> pmap;
        boost::tie (pmap, okay) = ps.template property_map<boost::uint8_t>(prop[i]);
        if (okay)
          {
            os << "property uchar " << prop[i] << std::endl;
            printers.push_back (new internal::Char_property_printer<Point,Vector,boost::uint8_t>(pmap));
            continue;
          }
      }
      {
        typename Point_set::template Property_map<boost::int16_t> pmap;
        boost::tie (pmap, okay) = ps.template property_map<boost::int16_t>(prop[i]);
        if (okay)
          {
            os << "property short " << prop[i] << std::endl;
            printers.push_back (new internal::Simple_property_printer<Point,Vector,boost::int16_t>(pmap));
            continue;
          }
      }
      {
        typename Point_set::template Property_map<boost::uint16_t> pmap;
        boost::tie (pmap, okay) = ps.template property_map<boost::uint16_t>(prop[i]);
        if (okay)
          {
            os << "property ushort " << prop[i] << std::endl;
            printers.push_back (new internal::Simple_property_printer<Point,Vector,boost::uint16_t>(pmap));
            continue;
          }
      }
      {
        typename Point_set::template Property_map<boost::int32_t> pmap;
        boost::tie (pmap, okay) = ps.template property_map<boost::int32_t>(prop[i]);
        if (okay)
          {
            os << "property int " << prop[i] << std::endl;
            printers.push_back (new internal::Simple_property_printer<Point,Vector,boost::int32_t>(pmap));
            continue;
          }
      }
      {
        typename Point_set::template Property_map<float> pmap;
        boost::tie (pmap, okay) = ps.template property_map<float>(prop[i]);
        if (okay)
          {
            os << "property float " << prop[i] << std::endl;
            printers.push_back (new internal::Simple_property_printer<Point,Vector,float>(pmap));
            continue;
          }
      }
      {
        typename Point_set::template Property_map<double> pmap;
        boost::tie (pmap, okay) = ps.template property_map<double>(prop[i]);
        if (okay)
          {
            os << "property double " << prop[i] << std::endl;
            printers.push_back (new internal::Simple_property_printer<Point,Vector,double>(pmap));
            continue;
          }
      }
    }
    
  os << "end_header" << std::endl;  

  for (typename Point_set::const_iterator it = ps.begin(); it != ps.end(); ++ it)
    {
      for (std::size_t i = 0; i < printers.size(); ++ i)
        {
          printers[i]->print(os, *it);
          if (get_mode (os) == IO::ASCII)
            os << " ";
        }
      if (get_mode (os) == IO::ASCII)
        os << std::endl;
    }
  return os;
}


  
} // namespace CGAL


#endif // CGAL_POINT_SET_3_IO
