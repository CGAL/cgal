// Copyright (c) 2016  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_3_IO
#define CGAL_POINT_SET_3_IO

#include <CGAL/license/Point_set_3.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/IO/write_off_points.h>

#include <CGAL/config.h>
#ifdef CGAL_LINKED_WITH_LASLIB
#include <CGAL/IO/read_las_points.h>
#include <CGAL/IO/write_las_points.h>
#endif // LAS
#include <CGAL/IO/PLY.h>

namespace CGAL {

namespace internal
{

namespace PLY
{

template <typename Point,
          typename Vector = typename Kernel_traits<Point>::Kernel::Vector_3>
class Point_set_3_filler
{
public:
  typedef Point_set_3<Point, Vector> Point_set;

private:

  struct Abstract_ply_property_to_point_set_property
  {
    virtual ~Abstract_ply_property_to_point_set_property() { }
    virtual void assign (PLY_element& element, typename Point_set::Index index) = 0;
  };

  template <typename Type>
  class PLY_property_to_point_set_property : public Abstract_ply_property_to_point_set_property
  {
    typedef typename Point_set::template Property_map<Type> Map;
    typedef typename Point_set::template Push_property_map<Map> Pmap;
    Map m_map;
    Pmap m_pmap;
    std::string m_name;
  public:
    PLY_property_to_point_set_property (Point_set& ps, const std::string& name)
      : m_name (name)
    {
      boost::tie (m_map, boost::tuples::ignore) = ps.add_property_map(name, Type());
      m_pmap = ps.push_property_map (m_map);
    }

    virtual void assign (PLY_element& element, typename Point_set::Index index)
    {
      Type t{};
      element.assign (t, m_name.c_str());
      put(m_pmap, index, t);
    }
  };

  Point_set& m_point_set;
  bool m_use_floats;
  std::vector<Abstract_ply_property_to_point_set_property*> m_properties;

public:

  Point_set_3_filler (Point_set& point_set)
    : m_point_set (point_set), m_use_floats (false)
  { }

  ~Point_set_3_filler()
  {
    for (std::size_t i = 0; i < m_properties.size(); ++ i)
      delete m_properties[i];
  }

  void instantiate_properties  (PLY_element& element)
  {
    bool has_normal[3] = { false, false, false };

    for (std::size_t j = 0; j < element.number_of_properties(); ++ j)
    {
      internal::PLY::PLY_read_number* property = element.property(j);

      const std::string& name = property->name();
      if (name == "x" ||
          name == "y" ||
          name == "z")
      {
        if (dynamic_cast<PLY_read_typed_number<float>*>(property))
          m_use_floats = true;
        continue;
      }
      if (name == "nx")
      {
        has_normal[0] = true;
        continue;
      }
      if (name == "ny")
      {
        has_normal[1] = true;
        continue;
      }
      if (name == "nz")
      {
        has_normal[2] = true;
        continue;
      }

      if (dynamic_cast<PLY_read_typed_number<boost::int8_t>*>(property))
      {
        m_properties.push_back
          (new PLY_property_to_point_set_property<boost::int8_t>(m_point_set,
                                                                 name));
      }
      else if (dynamic_cast<PLY_read_typed_number<boost::uint8_t>*>(property))
      {
        m_properties.push_back
          (new PLY_property_to_point_set_property<boost::uint8_t>(m_point_set,
                                                                  name));
      }
      else if (dynamic_cast<PLY_read_typed_number<boost::int16_t>*>(property))
      {
        m_properties.push_back
          (new PLY_property_to_point_set_property<boost::int16_t>(m_point_set,
                                                                  name));
      }
      else if (dynamic_cast<PLY_read_typed_number<boost::uint16_t>*>(property))
      {
        m_properties.push_back
          (new PLY_property_to_point_set_property<boost::uint16_t>(m_point_set,
                                                                   name));
      }
      else if (dynamic_cast<PLY_read_typed_number<boost::int32_t>*>(property))
      {
        m_properties.push_back
          (new PLY_property_to_point_set_property<boost::int32_t>(m_point_set,
                                                                  name));
      }
      else if (dynamic_cast<PLY_read_typed_number<boost::uint32_t>*>(property))
      {
        m_properties.push_back
          (new PLY_property_to_point_set_property<boost::uint32_t>(m_point_set,
                                                                   name));
      }
      else if (dynamic_cast<PLY_read_typed_number<float>*>(property))
      {
        m_properties.push_back
          (new PLY_property_to_point_set_property<float>(m_point_set,
                                                         name));
      }
      else if (dynamic_cast<PLY_read_typed_number<double>*>(property))
      {
        m_properties.push_back
          (new PLY_property_to_point_set_property<double>(m_point_set,
                                                          name));
      }
    }
    if (has_normal[0] && has_normal[1] && has_normal[2])
      m_point_set.add_normal_map();
  }

  void process_line (PLY_element& element)
  {
    m_point_set.insert();

    if (m_use_floats)
      process_line<float>(element);
    else
      process_line<double>(element);

    for (std::size_t i = 0; i < m_properties.size(); ++ i)
      m_properties[i]->assign (element, *(m_point_set.end() - 1));
  }

  template <typename FT>
  void process_line (PLY_element& element)
  {
    FT x = (FT)0.,y = (FT)0., z = (FT)0.,
      nx = (FT)0., ny = (FT)0., nz = (FT)0.;
    element.assign (x, "x");
    element.assign (y, "y");
    element.assign (z, "z");
    Point point (x, y, z);
    m_point_set.point(*(m_point_set.end() - 1)) = point;

    if (m_point_set.has_normal_map())
      {
        element.assign (nx, "nx");
        element.assign (ny, "ny");
        element.assign (nz, "nz");
        Vector normal (nx, ny, nz);
        m_point_set.normal(*(m_point_set.end() - 1)) = normal;
      }
  }
};

}

}

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

  bool out = CGAL::read_xyz_points
    (stream,
     point_set.index_back_inserter(),
     CGAL::parameters::point_map(point_set.point_push_map()).
     normal_map(point_set.normal_push_map()));

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

  bool out = CGAL::read_off_points
    (stream,
     point_set.index_back_inserter(),
     CGAL::parameters::point_map(point_set.point_push_map()).
     normal_map(point_set.normal_push_map()));

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


/// \cond SKIP_IN_MANUAL
template <typename Point, typename Vector>
bool
read_ply_point_set(
  std::istream& stream, ///< input stream.
  CGAL::Point_set_3<Point, Vector>& point_set) ///< point set
{
  std::string dummy;
  return read_ply_point_set (stream, point_set, dummy);
}

/// \endcond

/*!
  \ingroup PkgPointSet3IO

  Reads a point set with properties from an input stream in Ascii or
  Binary PLY format.

  - the operator reads the vertex `point` property;
  - if three PLY properties `nx`, `ny` and `nz` with type `float`
     or `double` are found, the normal map is added;
  - if any other PLY property is found, a "[name]" property map is
    added, where `[name]` is the name of the PLY property.

  The `comments` parameter can be omitted. If provided, it will be
  used to store the potential comments found in the PLY
  header. Each line starting by "comment " in the header is
  appended to the `comments` string (without the "comment " word).
 */
template <typename Point, typename Vector>
bool
read_ply_point_set(
  std::istream& stream, ///< input stream.
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  std::string& comments) ///< PLY comments.
{
  if(!stream)
    {
      std::cerr << "Error: cannot open file" << std::endl;
      return false;
    }

  internal::PLY::PLY_reader reader;
  internal::PLY::Point_set_3_filler<Point, Vector> filler(point_set);

  if (!(reader.init (stream)))
  {
    stream.setstate(std::ios::failbit);
    return false;
  }

  comments = reader.comments();

  for (std::size_t i = 0; i < reader.number_of_elements(); ++ i)
  {
    internal::PLY::PLY_element& element = reader.element(i);

    bool is_vertex = (element.name() == "vertex" || element.name() == "vertices");
    if (is_vertex)
    {
      point_set.reserve (element.number_of_items());
      filler.instantiate_properties (element);
    }

    for (std::size_t j = 0; j < element.number_of_items(); ++ j)
    {
      for (std::size_t k = 0; k < element.number_of_properties(); ++ k)
      {
        internal::PLY::PLY_read_number* property = element.property(k);
        property->get (stream);
        if (stream.fail())
          return false;
      }

      if (is_vertex)
        filler.process_line (element);
    }
  }

  return !stream.bad();
}

/*!
  \ingroup PkgPointSet3IO

  Writes a point set with properties in an output stream in PLY
  format.

  If found, the normal map is inserted to the stream.  All other
  properties with simple types are inserted in the stream.

  If provided, the `comments` string is included line by line in
  the header of the PLY stream (each line will be precedeed by
  "comment ").
 */
template <typename Point, typename Vector>
bool
write_ply_point_set(
  std::ostream& stream, ///< output stream.
  const CGAL::Point_set_3<Point, Vector>& point_set,  ///< point set.
  const std::string& comments = std::string()) ///< PLY comments.
{
  typedef CGAL::Point_set_3<Point, Vector> Point_set;
  typedef typename Point_set::Index Index;
  typedef typename Point_set::Point_map Point_map;
  typedef typename Point_set::Vector_map Vector_map;
  typedef typename Point_set::template Property_map<boost::int8_t> Int8_map;
  typedef typename Point_set::template Property_map<boost::uint8_t> Uint8_map;
  typedef typename Point_set::template Property_map<boost::int16_t> Int16_map;
  typedef typename Point_set::template Property_map<boost::uint16_t> Uint16_map;
  typedef typename Point_set::template Property_map<boost::int32_t> Int32_map;
  typedef typename Point_set::template Property_map<boost::uint32_t> Uint32_map;
  typedef typename Point_set::template Property_map<boost::int64_t> Int64_map;
  typedef typename Point_set::template Property_map<boost::uint64_t> Uint64_map;
  typedef typename Point_set::template Property_map<float> Float_map;
  typedef typename Point_set::template Property_map<double> Double_map;

  stream << "ply" << std::endl
         << ((get_mode(stream) == IO::BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
         << "comment Generated by the CGAL library" << std::endl;

  if (comments != std::string())
  {
    std::istringstream iss (comments);
    std::string line;
    while (getline(iss, line))
    {
      if (line != "Generated by the CGAL library") // Avoid repeating the line if multiple savings
        stream << "comment " << line << std::endl;
    }
  }

  stream << "element vertex " << point_set.number_of_points() << std::endl;

  std::vector<std::string> prop = point_set.base().properties();
  std::vector<internal::PLY::Abstract_property_printer<Index>*> printers;

  for (std::size_t i = 0; i < prop.size(); ++ i)
    {
      if (prop[i] == "index")
        continue;

      if (prop[i] == "point")
        {
          if (boost::is_same<typename Get_FT_from_map<typename Point_set::Point_map>::type, float>::value)
          {
            stream << "property float x" << std::endl
                   << "property float y" << std::endl
                   << "property float z" << std::endl;
          }
          else
          {
            stream << "property double x" << std::endl
                   << "property double y" << std::endl
                   << "property double z" << std::endl;
          }
          printers.push_back (new internal::PLY::Property_printer<Index,Point_map>(point_set.point_map()));
          continue;
        }
      if (prop[i] == "normal")
        {
          if (boost::is_same<typename Get_FT_from_map<typename Point_set::Vector_map>::type, float>::value)
          {
            stream << "property float nx" << std::endl
                   << "property float ny" << std::endl
                   << "property float nz" << std::endl;
          }
          else
          {
            stream << "property double nx" << std::endl
                   << "property double ny" << std::endl
                   << "property double nz" << std::endl;
          }
          printers.push_back (new internal::PLY::Property_printer<Index,Vector_map>(point_set.normal_map()));
          continue;
        }

      bool okay = false;
      {
        Int8_map pmap;
        boost::tie (pmap, okay) = point_set.template property_map<boost::int8_t>(prop[i]);
        if (okay)
          {
            stream << "property char " << prop[i] << std::endl;
            printers.push_back (new internal::PLY::Char_property_printer<Index,Int8_map>(pmap));
            continue;
          }
      }
      {
        Uint8_map pmap;
        boost::tie (pmap, okay) = point_set.template property_map<boost::uint8_t>(prop[i]);
        if (okay)
          {
            stream << "property uchar " << prop[i] << std::endl;
            printers.push_back (new internal::PLY::Char_property_printer<Index,Uint8_map>(pmap));
            continue;
          }
      }
      {
        Int16_map pmap;
        boost::tie (pmap, okay) = point_set.template property_map<boost::int16_t>(prop[i]);
        if (okay)
          {
            stream << "property short " << prop[i] << std::endl;
            printers.push_back (new internal::PLY::Simple_property_printer<Index,Int16_map>(pmap));
            continue;
          }
      }
      {
        Uint16_map pmap;
        boost::tie (pmap, okay) = point_set.template property_map<boost::uint16_t>(prop[i]);
        if (okay)
          {
            stream << "property ushort " << prop[i] << std::endl;
            printers.push_back (new internal::PLY::Simple_property_printer<Index,Uint16_map>(pmap));
            continue;
          }
      }
      {
        Int32_map pmap;
        boost::tie (pmap, okay) = point_set.template property_map<boost::int32_t>(prop[i]);
        if (okay)
          {
            stream << "property int " << prop[i] << std::endl;
            printers.push_back (new internal::PLY::Simple_property_printer<Index,Int32_map>(pmap));
            continue;
          }
      }
      {
        Uint32_map pmap;
        boost::tie (pmap, okay) = point_set.template property_map<boost::uint32_t>(prop[i]);
        if (okay)
          {
            stream << "property uint " << prop[i] << std::endl;
            printers.push_back (new internal::PLY::Simple_property_printer<Index,Uint32_map>(pmap));
            continue;
          }
      }
      {
        Int64_map pmap;
        boost::tie (pmap, okay) = point_set.template property_map<boost::int64_t>(prop[i]);
        if (okay)
          {
            stream << "property int " << prop[i] << std::endl;
            printers.push_back (new internal::PLY::Simple_property_printer<Index,Int64_map,boost::int32_t>(pmap));
            continue;
          }
      }
      {
        Uint64_map pmap;
        boost::tie (pmap, okay) = point_set.template property_map<boost::uint64_t>(prop[i]);
        if (okay)
          {
            stream << "property uint " << prop[i] << std::endl;
            printers.push_back (new internal::PLY::Simple_property_printer<Index,Uint64_map,boost::uint32_t>(pmap));
            continue;
          }
      }
      {
        Float_map pmap;
        boost::tie (pmap, okay) = point_set.template property_map<float>(prop[i]);
        if (okay)
          {
            stream << "property float " << prop[i] << std::endl;
            printers.push_back (new internal::PLY::Simple_property_printer<Index,Float_map>(pmap));
            continue;
          }
      }
      {
        Double_map pmap;
        boost::tie (pmap, okay) = point_set.template property_map<double>(prop[i]);
        if (okay)
          {
            stream << "property double " << prop[i] << std::endl;
            printers.push_back (new internal::PLY::Simple_property_printer<Index,Double_map>(pmap));
            continue;
          }
      }
    }

  stream << "end_header" << std::endl;

  for (typename Point_set::const_iterator it = point_set.begin(); it != point_set.end(); ++ it)
    {
      for (std::size_t i = 0; i < printers.size(); ++ i)
        {
          printers[i]->print(stream, *it);
          if (get_mode (stream) == IO::ASCII)
            stream << " ";
        }
      if (get_mode (stream) == IO::ASCII)
        stream << std::endl;
    }

  for (std::size_t i = 0; i < printers.size(); ++ i)
    delete printers[i];
  return true;
}

#if defined(CGAL_LINKED_WITH_LASLIB) || defined(DOXYGEN_RUNNING)

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

/*!
  \ingroup PkgPointSet3IO
 */
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
    = read_las_points_with_properties
    (stream, point_set.index_back_inserter(),
     make_las_point_reader (point_set.point_push_map()),
     std::make_pair (point_set.push_property_map (intensity), LAS_property::Intensity()),
     std::make_pair (point_set.push_property_map (return_number), LAS_property::Return_number()),
     std::make_pair (point_set.push_property_map (number_of_returns), LAS_property::Number_of_returns()),
     std::make_pair (point_set.push_property_map (scan_direction_flag), LAS_property::Scan_direction_flag()),
     std::make_pair (point_set.push_property_map (edge_of_flight_line), LAS_property::Edge_of_flight_line()),
     std::make_pair (point_set.push_property_map (classification), LAS_property::Classification()),
     std::make_pair (point_set.push_property_map (synthetic_flag), LAS_property::Synthetic_flag()),
     std::make_pair (point_set.push_property_map (keypoint_flag), LAS_property::Keypoint_flag()),
     std::make_pair (point_set.push_property_map (withheld_flag), LAS_property::Withheld_flag()),
     std::make_pair (point_set.push_property_map (scan_angle), LAS_property::Scan_angle()),
     std::make_pair (point_set.push_property_map (user_data), LAS_property::User_data()),
     std::make_pair (point_set.push_property_map (point_source_ID), LAS_property::Point_source_ID()),
     std::make_pair (point_set.push_property_map (deleted_flag), LAS_property::Deleted_flag()),
     std::make_pair (point_set.push_property_map (gps_time), LAS_property::GPS_time()),
     std::make_pair (point_set.push_property_map (R), LAS_property::R()),
     std::make_pair (point_set.push_property_map (G), LAS_property::G()),
     std::make_pair (point_set.push_property_map (B), LAS_property::B()),
     std::make_pair (point_set.push_property_map (I), LAS_property::I()));

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
  internal::check_if_property_is_used (point_set, gps_time);
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

  Double_map gps_time;
  bool remove_gps_time;
  boost::tie (gps_time, remove_gps_time)
    = point_set.template add_property_map<double>("gps_time", 0);

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
    (stream, point_set,
     make_las_point_writer (point_set.point_map()),
     std::make_pair (intensity, LAS_property::Intensity()),
     std::make_pair (return_number, LAS_property::Return_number()),
     std::make_pair (number_of_returns, LAS_property::Number_of_returns()),
     std::make_pair (scan_direction_flag, LAS_property::Scan_direction_flag()),
     std::make_pair (edge_of_flight_line, LAS_property::Edge_of_flight_line()),
     std::make_pair (classification, LAS_property::Classification()),
     std::make_pair (synthetic_flag, LAS_property::Synthetic_flag()),
     std::make_pair (keypoint_flag, LAS_property::Keypoint_flag()),
     std::make_pair (withheld_flag, LAS_property::Withheld_flag()),
     std::make_pair (scan_angle, LAS_property::Scan_angle()),
     std::make_pair (user_data, LAS_property::User_data()),
     std::make_pair (point_source_ID, LAS_property::Point_source_ID()),
     std::make_pair (deleted_flag, LAS_property::Deleted_flag()),
     std::make_pair (gps_time, LAS_property::GPS_time()),
     std::make_pair (R, LAS_property::R()),
     std::make_pair (G, LAS_property::G()),
     std::make_pair (B, LAS_property::B()),
     std::make_pair (I, LAS_property::I()));

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
  if (remove_gps_time) point_set.remove_property_map (gps_time);
  if (remove_R) point_set.remove_property_map (R);
  if (remove_G) point_set.remove_property_map (G);
  if (remove_B) point_set.remove_property_map (B);
  if (remove_I) point_set.remove_property_map (I);

  return okay;
}

#endif // LAS

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
    return CGAL::write_xyz_points
      (stream, point_set,
       CGAL::parameters::point_map(point_set.point_map()).
       normal_map(point_set.normal_map()));

  return CGAL::write_xyz_points
    (stream, point_set,
     CGAL::parameters::point_map(point_set.point_map()));
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
    return CGAL::write_off_points
      (stream, point_set,
       CGAL::parameters::point_map(point_set.point_map()).
       normal_map(point_set.normal_map()));

  return CGAL::write_off_points
    (stream, point_set,
     CGAL::parameters::point_map(point_set.point_map()));
}


/*!

  \ingroup PkgPointSet3IO

  \brief Reads the point set from an input stream that can be either:

  - XYZ
  - OFF
  - PLY
  - LAS

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
  else if (line.find("LASF") == 0)
    CGAL::read_las_point_set (is, ps);
#endif // LAS
  else
    CGAL::read_xyz_point_set (is, ps);

  return is;
}

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
  write_ply_point_set (os, ps);
  return os;
}



} // namespace CGAL


#endif // CGAL_POINT_SET_3_IO
