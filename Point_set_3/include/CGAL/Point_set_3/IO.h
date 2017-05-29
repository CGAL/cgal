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
  // if (point_set.has_normal_map())
  //   return CGAL::write_ply_points_and_normals
  //     (stream, point_set.begin(), point_set.end(),
  //      point_set.point_map(), point_set.normal_map());
  
  // return CGAL::write_ply_points
  // (stream, point_set.begin(), point_set.end(),
  //  point_set.point_map());
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
    virtual std::string get_string(const typename CGAL::Point_set_3<Point,Vector>::Index& index) = 0;
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
    
    virtual std::string get_string(const typename CGAL::Point_set_3<Point,Vector>::Index& index)
    {
      std::ostringstream oss;
      oss.precision (std::numeric_limits<double>::digits10 + 2);
      oss << get(m_pmap, index);
      return oss.str();
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
    
    virtual std::string get_string(const typename CGAL::Point_set_3<Point,Vector>::Index& index)
    {
      std::ostringstream oss;
      oss << (int)(get(m_pmap, index));
      return oss.str();
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
     << "format ascii 1.0" << std::endl
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
            printers.push_back (new internal::Property_printer<Point,Vector,boost::int16_t>(pmap));
            continue;
          }
      }
      {
        typename Point_set::template Property_map<boost::uint16_t> pmap;
        boost::tie (pmap, okay) = ps.template property_map<boost::uint16_t>(prop[i]);
        if (okay)
          {
            os << "property ushort " << prop[i] << std::endl;
            printers.push_back (new internal::Property_printer<Point,Vector,boost::uint16_t>(pmap));
            continue;
          }
      }
      {
        typename Point_set::template Property_map<boost::int32_t> pmap;
        boost::tie (pmap, okay) = ps.template property_map<boost::int32_t>(prop[i]);
        if (okay)
          {
            os << "property int " << prop[i] << std::endl;
            printers.push_back (new internal::Property_printer<Point,Vector,boost::int32_t>(pmap));
            continue;
          }
      }
      {
        typename Point_set::template Property_map<float> pmap;
        boost::tie (pmap, okay) = ps.template property_map<float>(prop[i]);
        if (okay)
          {
            os << "property float " << prop[i] << std::endl;
            printers.push_back (new internal::Property_printer<Point,Vector,float>(pmap));
            continue;
          }
      }
      {
        typename Point_set::template Property_map<double> pmap;
        boost::tie (pmap, okay) = ps.template property_map<double>(prop[i]);
        if (okay)
          {
            os << "property double " << prop[i] << std::endl;
            printers.push_back (new internal::Property_printer<Point,Vector,double>(pmap));
            continue;
          }
      }
    }
    
  os << "end_header" << std::endl;  

  for (typename Point_set::const_iterator it = ps.begin(); it != ps.end(); ++ it)
    {
      for (std::size_t i = 0; i < printers.size(); ++ i)
        os << printers[i]->get_string(*it) << " ";
      os << std::endl;
    }
  return os;
}


  
} // namespace CGAL


#endif // CGAL_POINT_SET_3_IO
