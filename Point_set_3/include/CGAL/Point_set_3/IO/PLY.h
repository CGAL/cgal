// Copyright (c) 2016 GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_IO_PLY_H
#define CGAL_POINT_SET_IO_PLY_H

#include <CGAL/license/Point_set_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/IO/io.h>

#include <fstream>
#include <string>
#include <vector>

namespace CGAL {

template <typename Point, typename Vector>
class Point_set_3;

namespace IO {
namespace internal {

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
    virtual void assign(PLY_element& element, typename Point_set::Index index) = 0;
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
    PLY_property_to_point_set_property(Point_set& ps, const std::string& name)
      : m_name(name)
    {
      boost::tie(m_map, boost::tuples::ignore) = ps.add_property_map(name, Type());
      m_pmap = ps.push_property_map(m_map);
    }

    virtual void assign(PLY_element& element, typename Point_set::Index index)
    {
      Type t{};
      element.assign(t, m_name.c_str());
      put(m_pmap, index, t);
    }
  };

  Point_set& m_point_set;
  bool m_use_floats;
  std::vector<Abstract_ply_property_to_point_set_property*> m_properties;

public:

  Point_set_3_filler(Point_set& point_set)
    : m_point_set(point_set), m_use_floats(false)
  { }

  ~Point_set_3_filler()
  {
    for(std::size_t i=0; i<m_properties.size(); ++i)
      delete m_properties[i];
  }

  void instantiate_properties(PLY_element& element)
  {
    bool has_normal[3] = { false, false, false };

    for(std::size_t j=0; j<element.number_of_properties(); ++j)
    {
      internal::PLY_read_number* property = element.property(j);

      const std::string& name = property->name();
      if(name == "x" ||
          name == "y" ||
          name == "z")
      {
        if(dynamic_cast<PLY_read_typed_number<float>*>(property))
          m_use_floats = true;
        continue;
      }
      if(name == "nx")
      {
        has_normal[0] = true;
        continue;
      }
      if(name == "ny")
      {
        has_normal[1] = true;
        continue;
      }
      if(name == "nz")
      {
        has_normal[2] = true;
        continue;
      }

      if(dynamic_cast<PLY_read_typed_number<std::int8_t>*>(property))
      {
        m_properties.push_back
            (new PLY_property_to_point_set_property<std::int8_t>(m_point_set,
                                                                   name));
      }
      else if(dynamic_cast<PLY_read_typed_number<std::uint8_t>*>(property))
      {
        m_properties.push_back
            (new PLY_property_to_point_set_property<std::uint8_t>(m_point_set,
                                                                    name));
      }
      else if(dynamic_cast<PLY_read_typed_number<std::int16_t>*>(property))
      {
        m_properties.push_back
            (new PLY_property_to_point_set_property<std::int16_t>(m_point_set,
                                                                    name));
      }
      else if(dynamic_cast<PLY_read_typed_number<std::uint16_t>*>(property))
      {
        m_properties.push_back
            (new PLY_property_to_point_set_property<std::uint16_t>(m_point_set,
                                                                     name));
      }
      else if(dynamic_cast<PLY_read_typed_number<std::int32_t>*>(property))
      {
        m_properties.push_back
            (new PLY_property_to_point_set_property<std::int32_t>(m_point_set,
                                                                    name));
      }
      else if(dynamic_cast<PLY_read_typed_number<std::uint32_t>*>(property))
      {
        m_properties.push_back
            (new PLY_property_to_point_set_property<std::uint32_t>(m_point_set,
                                                                     name));
      }
      else if(dynamic_cast<PLY_read_typed_number<float>*>(property))
      {
        m_properties.push_back
            (new PLY_property_to_point_set_property<float>(m_point_set,
                                                           name));
      }
      else if(dynamic_cast<PLY_read_typed_number<double>*>(property))
      {
        m_properties.push_back
            (new PLY_property_to_point_set_property<double>(m_point_set,
                                                            name));
      }
    }
    if(has_normal[0] && has_normal[1] && has_normal[2])
      m_point_set.add_normal_map();
  }

  void process_line(PLY_element& element)
  {
    m_point_set.insert();

    if(m_use_floats)
      process_line<float>(element);
    else
      process_line<double>(element);

    for(std::size_t i=0; i<m_properties.size(); ++i)
      m_properties[i]->assign(element, *(m_point_set.end() - 1));
  }

  template <typename FT>
  void process_line(PLY_element& element)
  {
    FT x = (FT)0.,y = (FT)0., z = (FT)0.,
        nx = (FT)0., ny = (FT)0., nz = (FT)0.;
    element.assign(x, "x");
    element.assign(y, "y");
    element.assign(z, "z");
    Point point(x, y, z);
    m_point_set.point(*(m_point_set.end() - 1)) = point;

    if(m_point_set.has_normal_map())
    {
      element.assign(nx, "nx");
      element.assign(ny, "ny");
      element.assign(nz, "nz");
      Vector normal(nx, ny, nz);
      m_point_set.normal(*(m_point_set.end() - 1)) = normal;
    }
  }
};

} // namespace internal

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

/*!
  \ingroup PkgPointSet3IOPLY

  \brief reads a point set with properties from an input stream in \ascii or binary \ref IOStreamPLY.

  - the operator reads the vertex `point` property;
  - if three PLY properties `nx`, `ny` and `nz` with type `float`
     or `double` are found, the normal map is added;
  - if any other PLY property is found, a "[name]" property map is
    added, where `[name]` is the name of the PLY property.

  The `comments` parameter can be omitted. If provided, it will be
  used to store the potential comments found in the PLY
  header. Each line starting by "comment " in the header is
  appended to the `comments` string (without the "comment " word).

  \attention To read a binary file, the flag `std::ios::binary` must be set during the creation of the `ifstream`.

  \param is the input stream
  \param point_set the point set
  \param comments optional PLY comments.

  \return `true` if the reading was successful, `false` otherwise.
 */
template <typename Point, typename Vector>
bool read_PLY(std::istream& is,
              CGAL::Point_set_3<Point, Vector>& point_set,
              std::string& comments)
{
  if(!is)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  internal::PLY_reader reader(true);
  internal::Point_set_3_filler<Point, Vector> filler(point_set);

  if(!(reader.init(is)))
  {
    is.setstate(std::ios::failbit);
    return false;
  }

  comments = reader.comments();

  for(std::size_t i=0; i<reader.number_of_elements(); ++i)
  {
    internal::PLY_element& element = reader.element(i);

    bool is_vertex = (element.name() == "vertex" || element.name() == "vertices");
    if(is_vertex)
    {
      point_set.reserve(element.number_of_items());
      filler.instantiate_properties(element);
    }

    for(std::size_t j=0; j<element.number_of_items(); ++j)
    {
      for(std::size_t k=0; k<element.number_of_properties(); ++k)
      {
        internal::PLY_read_number* property = element.property(k);
        property->get(is);
        if(is.fail())
          return false;
      }

      if(is_vertex)
        filler.process_line(element);
    }
  }

  return !is.bad();
}

/// \cond SKIP_IN_MANUAL

template <typename Point, typename Vector>
bool read_PLY(std::istream& is, CGAL::Point_set_3<Point, Vector>& point_set)
{
  std::string dummy;
  return read_PLY(is, point_set, dummy);
}

/// \endcond

/*!
  \ingroup PkgPointSet3IOPLY

  \brief reads a point set with properties from an input stream in \ascii or binary \ref IOStreamPLY.

  - the operator reads the vertex `point` property;
  - if three PLY properties `nx`, `ny` and `nz` with type `float`
     or `double` are found, the normal map is added;
  - if any other PLY property is found, a "[name]" property map is
    added, where `[name]` is the name of the PLY property.

  The `comments` parameter can be omitted. If provided, it will be
  used to store the potential comments found in the PLY
  header. Each line starting by "comment " in the header is
  appended to the `comments` string (without the "comment " word).

  \tparam Point the point type of the `Point_set_3`
  \tparam Vector the vector type of the `Point_set_3`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the path to the input file
  \param point_set the point set
  \param comments optional PLY comments.
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamNBegin{use_binary_mode}
      \cgalParamDescription{indicates whether data should be read in binary (`true`) or in \ascii (`false`)}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \return `true` if the reading was successful, `false` otherwise.
*/
template <typename Point, typename Vector, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const std::string& fname,
              CGAL::Point_set_3<Point, Vector>& point_set,
              std::string& comments,
              const CGAL_NP_CLASS& np = parameters::default_values())
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ifstream is(fname, std::ios::binary);
    CGAL::IO::set_mode(is, CGAL::IO::BINARY);
    return read_PLY(is, point_set, comments);
  }
  else
  {
    std::ifstream is(fname);
    CGAL::IO::set_mode(is, CGAL::IO::ASCII);
    return read_PLY(is, point_set, comments);
  }
}

/// \cond SKIP_IN_MANUAL
template <typename Point, typename Vector, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const std::string& fname, CGAL::Point_set_3<Point, Vector>& point_set, const CGAL_NP_CLASS& np = parameters::default_values())
{
  std::string unused_comments;
  return read_PLY(fname, point_set, unused_comments, np);
}
/// \endcond

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgPointSet3IODeprecated

  \deprecated This function is deprecated since \cgal 5.3,
              \link PkgPointSet3IO `CGAL::IO::read_PLY()` \endlink  should be used instead.

  \brief reads a point set with properties from an input stream in \ascii or binary PLY format.

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
CGAL_DEPRECATED bool read_ply_point_set(std::istream& is, ///< input stream.
                                        CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
                                        std::string& comments) ///< PLY comments.
{
  return IO::read_PLY(is, point_set, comments);
}

template <typename Point, typename Vector>
CGAL_DEPRECATED bool read_ply_point_set(std::istream& is, ///< input stream.
                                        CGAL::Point_set_3<Point, Vector>& point_set) ///< point set
{
  std::string dummy;
  return IO::read_PLY(is, point_set, dummy);
}

#endif // CGAL_NO_DEPRECATED_CODE

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace IO {

/*!
  \ingroup PkgPointSet3IOPLY

  \brief writes a point set with properties in an output stream in the \ref IOStreamPLY.

  If it exists, the normal map is inserted in the stream. All other
  properties with simple types are inserted in the stream.

  If provided, the `comments` string is included line by line in
  the header of the PLY stream (each line will be precedeed by
  "comment ").

  \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation
             of the `ofstream`, and the \link PkgStreamSupportEnumRef `IO::Mode` \endlink
             of the stream must be set to `BINARY`.

  \tparam Point the point type of the `Point_set_3`
  \tparam Vector the vector type of the `Point_set_3`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param os the output stream
  \param point_set the point set
  \param comments optional PLY comments
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
   \cgalParamNBegin{stream_precision}
     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
     \cgalParamType{int}
     \cgalParamDefault{the precision of the stream `os`}
     \cgalParamExtra{This parameter is only meaningful while using \ascii encoding.}
   \cgalParamNEnd
  \cgalNamedParamsEnd

  \return `true` if the reading was successful, `false` otherwise.
*/
template <typename Point, typename Vector, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(std::ostream& os,
               const CGAL::Point_set_3<Point, Vector>& point_set,
               const std::string& comments,
               const CGAL_NP_CLASS& np = parameters::default_values())
{
  typedef CGAL::Point_set_3<Point, Vector> Point_set;
  typedef typename Point_set::Index Index;
  typedef typename Point_set::Point_map Point_map;
  typedef typename Point_set::Vector_map Vector_map;
  typedef typename Point_set::template Property_map<std::int8_t> Int8_map;
  typedef typename Point_set::template Property_map<std::uint8_t> Uint8_map;
  typedef typename Point_set::template Property_map<std::int16_t> Int16_map;
  typedef typename Point_set::template Property_map<std::uint16_t> Uint16_map;
  typedef typename Point_set::template Property_map<std::int32_t> Int32_map;
  typedef typename Point_set::template Property_map<std::uint32_t> Uint32_map;
  typedef typename Point_set::template Property_map<std::int64_t> Int64_map;
  typedef typename Point_set::template Property_map<std::uint64_t> Uint64_map;
  typedef typename Point_set::template Property_map<float> Float_map;
  typedef typename Point_set::template Property_map<double> Double_map;

  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  if(!os)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  set_stream_precision_from_NP(os, np);

  os << "ply" << std::endl
     << ((CGAL::IO::get_mode(os) == CGAL::IO::BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
     << "comment Generated by the CGAL library" << std::endl;

  if(comments != std::string())
  {
    std::istringstream iss(comments);
    std::string line;
    while(getline(iss, line))
    {
      if(line != "Generated by the CGAL library") // Avoid repeating the line if multiple savings
        os << "comment " << line << std::endl;
    }
  }

  os << "element vertex " << point_set.number_of_points() << std::endl;

  std::vector<std::string> prop = point_set.base().properties();
  std::vector<internal::Abstract_property_printer<Index>*> printers;

  for(std::size_t i=0; i<prop.size(); ++i)
  {
    if(prop[i] == "index")
      continue;

    if(prop[i] == "point")
    {
      if(std::is_same<typename Get_FT_from_map<typename Point_set::Point_map>::type, float>::value)
      {
        os << "property float x" << std::endl
           << "property float y" << std::endl
           << "property float z" << std::endl;
      }
      else
      {
        os << "property double x" << std::endl
           << "property double y" << std::endl
           << "property double z" << std::endl;
      }
      printers.push_back(new internal::Property_printer<Index,Point_map>(point_set.point_map()));
      continue;
    }
    if(prop[i] == "normal")
    {
      if(std::is_same<typename Get_FT_from_map<typename Point_set::Vector_map>::type, float>::value)
      {
        os << "property float nx" << std::endl
           << "property float ny" << std::endl
           << "property float nz" << std::endl;
      }
      else
      {
        os << "property double nx" << std::endl
           << "property double ny" << std::endl
           << "property double nz" << std::endl;
      }
      printers.push_back(new internal::Property_printer<Index,Vector_map>(point_set.normal_map()));
      continue;
    }

    bool okay = false;
    {
      Int8_map pmap;
      boost::tie(pmap, okay) = point_set.template property_map<std::int8_t>(prop[i]);
      if(okay)
      {
        os << "property char " << prop[i] << std::endl;
        printers.push_back(new internal::Simple_property_printer<Index,Int8_map>(pmap));
        continue;
      }
    }
    {
      Uint8_map pmap;
      boost::tie(pmap, okay) = point_set.template property_map<std::uint8_t>(prop[i]);
      if(okay)
      {
        os << "property uchar " << prop[i] << std::endl;
        printers.push_back(new internal::Simple_property_printer<Index,Uint8_map>(pmap));
        continue;
      }
    }
    {
      Int16_map pmap;
      boost::tie(pmap, okay) = point_set.template property_map<std::int16_t>(prop[i]);
      if(okay)
      {
        os << "property short " << prop[i] << std::endl;
        printers.push_back(new internal::Simple_property_printer<Index,Int16_map>(pmap));
        continue;
      }
    }
    {
      Uint16_map pmap;
      boost::tie(pmap, okay) = point_set.template property_map<std::uint16_t>(prop[i]);
      if(okay)
      {
        os << "property ushort " << prop[i] << std::endl;
        printers.push_back(new internal::Simple_property_printer<Index,Uint16_map>(pmap));
        continue;
      }
    }
    {
      Int32_map pmap;
      boost::tie(pmap, okay) = point_set.template property_map<std::int32_t>(prop[i]);
      if(okay)
      {
        os << "property int " << prop[i] << std::endl;
        printers.push_back(new internal::Simple_property_printer<Index,Int32_map>(pmap));
        continue;
      }
    }
    {
      Uint32_map pmap;
      boost::tie(pmap, okay) = point_set.template property_map<std::uint32_t>(prop[i]);
      if(okay)
      {
        os << "property uint " << prop[i] << std::endl;
        printers.push_back(new internal::Simple_property_printer<Index,Uint32_map>(pmap));
        continue;
      }
    }
    {
      Int64_map pmap;
      boost::tie(pmap, okay) = point_set.template property_map<std::int64_t>(prop[i]);
      if(okay)
      {
        os << "property int " << prop[i] << std::endl;
        printers.push_back(new internal::Simple_property_printer<Index,Int64_map,std::int32_t>(pmap));
        continue;
      }
    }
    {
      Uint64_map pmap;
      boost::tie(pmap, okay) = point_set.template property_map<std::uint64_t>(prop[i]);
      if(okay)
      {
        os << "property uint " << prop[i] << std::endl;
        printers.push_back(new internal::Simple_property_printer<Index,Uint64_map,std::uint32_t>(pmap));
        continue;
      }
    }
    {
      Float_map pmap;
      boost::tie(pmap, okay) = point_set.template property_map<float>(prop[i]);
      if(okay)
      {
        os << "property float " << prop[i] << std::endl;
        printers.push_back(new internal::Simple_property_printer<Index,Float_map>(pmap));
        continue;
      }
    }
    {
      Double_map pmap;
      boost::tie(pmap, okay) = point_set.template property_map<double>(prop[i]);
      if(okay)
      {
        os << "property double " << prop[i] << std::endl;
        printers.push_back(new internal::Simple_property_printer<Index,Double_map>(pmap));
        continue;
      }
    }
  }

  os << "end_header" << std::endl;

  for(typename Point_set::const_iterator it = point_set.begin(); it != point_set.end(); ++it)
  {
    for(std::size_t i=0; i<printers.size(); ++i)
    {
      printers[i]->print(os, *it);
      if(get_mode(os) == ASCII)
        os << " ";
    }

    if(get_mode(os) == ASCII)
      os << std::endl;
  }

  for(std::size_t i=0; i<printers.size(); ++i)
    delete printers[i];
  return true;
}

/// \cond SKIP_IN_MANUAL

template <typename Point, typename Vector, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(std::ostream& os, const CGAL::Point_set_3<Point, Vector>& point_set, const CGAL_NP_CLASS& np = parameters::default_values())
{
  return write_PLY(os, point_set, std::string(), np);
}

/// \endcond

/*!
  \ingroup PkgPointSet3IOPLY

  \brief writes a point set with properties in an output stream in the \ref IOStreamPLY.

  If it exists, the normal map is written in the file. All other
  properties with simple types are written in the file.

  If provided, the `comments` string is included line by line in
  the header of the PLY stream (each line will be precedeed by
  "comment ").

  \tparam Point the point type of the `Point_set_3`
  \tparam Vector the vector type of the `Point_set_3`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the path to the output file
  \param point_set the point set
  \param comments optional PLY comments
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{use_binary_mode}
      \cgalParamDescription{indicates whether data should be read in binary (`true`) or in \ascii (`false`)}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
    \cgalParamNEnd

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`6`}
      \cgalParamExtra{This parameter is only meaningful while using \ascii encoding.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \return `true` if the reading was successful, `false` otherwise.
*/
template <typename Point, typename Vector, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(const std::string& fname,
               const CGAL::Point_set_3<Point, Vector>& point_set,
               const std::string& comments,
               const CGAL_NP_CLASS& np)
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ofstream os(fname, std::ios::binary);
    CGAL::IO::set_mode(os, BINARY);
    return write_PLY(os, point_set, comments, np);
  }
  else
  {
    std::ofstream os(fname);
    CGAL::IO::set_mode(os, ASCII);
    return write_PLY(os, point_set, comments, np);
  }
}

/// \cond SKIP_IN_MANUAL

template <typename Point, typename Vector, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(const std::string& fname, const CGAL::Point_set_3<Point, Vector>& point_set, const CGAL_NP_CLASS& np = parameters::default_values())
{
  return write_PLY(fname, point_set, std::string(), np);
}

/// \endcond

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgPointSet3IODeprecated

  \deprecated This function is deprecated since \cgal 5.3,
              \link PkgPointSet3IO `CGAL::IO::write_PLY()` \endlink  should be used instead.
 */
template <typename Point, typename Vector>
CGAL_DEPRECATED bool write_ply_point_set(std::ostream& os,
                                         const CGAL::Point_set_3<Point, Vector>& point_set,
                                         const std::string& comments = std::string())
{
  return IO::write_PLY(os, point_set, comments);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_POINT_SET_IO_PLY_H
