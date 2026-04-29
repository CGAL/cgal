// Copyright (c) 2020 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//                 Mael Rouxel-Labbé
//                 Maxime Gimeno

#ifndef CGAL_IO_GOCAD_H
#define CGAL_IO_GOCAD_H

#include <CGAL/IO/helpers.h>
#include <CGAL/IO/io.h>

#include <CGAL/assertions.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/use.h>

#include <boost/range/value_type.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <type_traits>

namespace CGAL {

namespace IO {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(std::istream& is,
                std::pair<std::string, std::string>& name_and_color,
                PointRange& points,
                PolygonRange& polygons,
                const CGAL_NP_CLASS& np)
{
  typedef typename boost::range_value<PointRange>::type     Point;
  typedef typename boost::range_value<PolygonRange>::type   Poly;

  const bool verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);
  const bool read_only_one_object = parameters::choose_parameter(parameters::get_parameter(np, internal_np::read_only_one_object), false);


  if(!is)
  {
    if(verbose)
      std::cerr<<"File doesn't exist."<<std::endl;
    return false;
  }

  set_ascii_mode(is); // GOCAD is ASCII only

  std::string s;
  Point p;

  bool vertices_read = false;
  std::size_t offset=0;
  int l_offset=0;
  int nb_gocad=0;
  int nb_end=0;
  bool end_read=true;

  std::string line;
  while(std::getline(is, line))
  {
    if(line.empty())
      continue;

    std::istringstream iss(line);
    if(!(iss >> s))
      continue; // can't read anything on the line, whitespace only?

    if(s == "GOCAD")
    {
      ++nb_gocad;
      offset=points.size();
      vertices_read = false;
      end_read = false;
    }
    else if(s=="VRTX")
    {
      int i;
      double x, y, z;
      if(!(iss >> i >> IO::iformat(x) >> IO::iformat(y) >> IO::iformat(z)))
      {
        if(verbose)
          std::cerr << "error while reading vertex." << std::endl;
        return false;
      }

      if(!vertices_read)
      {
        vertices_read = true;
        l_offset = -i;
      }

      internal::fill_point(x, y, z, 1., p);
      points.push_back(p);
    }
    else if(s == "TRGL")
    {
      int i,j,k;
      if(!(iss >> i >> j >> k))
      {
        if(verbose)
          std::cerr << "error while reading triangle." << std::endl;
        return false;
      }

      Poly new_face;
      ::CGAL::internal::resize(new_face, 3);
      new_face[0] =  offset + (l_offset + i);
      new_face[1] =  offset + (l_offset + j);
      new_face[2] =  offset + (l_offset + k);
      polygons.push_back(new_face);
    }
    else if (s=="END")
    {
      ++nb_end;
      end_read=true;
      if (read_only_one_object)
        break;
    }
    else if (s=="HEADER")
    {
      while(std::getline(is, line))
      {
        std::string::size_type idx;
        if((idx = line.find("name")) != std::string::npos)
        {
          std::size_t pos = line.find(":")+1;
          name_and_color.first = line.substr(pos, line.length());
        }
        else if((idx = line.find("color")) != std::string::npos)
        {
          std::size_t pos = line.find(":")+1;
          name_and_color.second = line.substr(pos, line.length());
        }
        else if((idx = line.find("}")) != std::string::npos)
          break;
      }
    }
  }

  return nb_gocad!=0 && end_read && nb_gocad == nb_end && (is.eof() || !is.bad());
}
/// \endcond

/*!
 * \ingroup PkgStreamSupportIoFuncsGOCAD
 *
 * \brief reads the content of `is` into `points` and `polygons`, using the \ref IOStreamGocad.
 *
 * \attention The polygon soup is not cleared, and the data from the stream are appended.
 *
 * \tparam PointRange a model of the concept `BackInsertionSequence` whose value type is the point type
 * \tparam PolygonRange a model of the concept `BackInsertionSequence` whose `value_type` is
 *                      a model of the concepts `RandomAccessContainer` whose `value_type` is
 *                      an unsigned integer type convertible to `std::size_t`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param is the input stream
 * \param points points of the soup of polygons
 * \param polygons a range of polygons. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{indicates whether output warnings and error messages should be printed or not.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{read_only_one_object}
 *     \cgalParamDescription{if `true` only one GOCAD object will be read in the stream,  and the stream will be left ready for the reading of the next GOCAD object.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the reading was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(std::istream& is,
                PointRange& points,
                PolygonRange& polygons,
                const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
                )
{
  std::pair<std::string, std::string> dummy;
  return read_GOCAD(is, dummy, points, polygons, np);
}

/*!
 * \ingroup PkgStreamSupportIoFuncsGOCAD
 *
 * \brief reads the content of the file `fname` into `points` and `polygons`, using the \ref IOStreamGocad.
 *
 * \attention The polygon soup is not cleared, and the data from the file are appended.
 *
 * \tparam PointRange a model of the concept `BackInsertionSequence` whose value type is the point type
 * \tparam PolygonRange a model of the concept `BackInsertionSequence` whose `value_type` is
 *                      a model of the concepts `RandomAccessContainer` whose `value_type` is
 *                      an unsigned integer type convertible to `std::size_t`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the path to the input file
 * \param points points of the soup of polygons
 * \param polygons a range of polygons. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{indicates whether output warnings and error messages should be printed or not.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{read_only_one_object}
 *     \cgalParamDescription{if `true` only one GOCAD object will be read in the stream,  and the stream will be left ready for the reading of the next GOCAD object.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the reading was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(const std::string& fname,
                PointRange& points,
                PolygonRange& polygons,
                const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
                )
{
  std::ifstream is(fname);
  CGAL::IO::set_mode(is, CGAL::IO::ASCII);
  std::pair<std::string, std::string> dummy;
  return read_GOCAD(is, dummy, points, polygons, np);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace internal {

template <typename PointRange,
          typename PolygonRange,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(std::ostream& os,
                 const char* fname,
                 const PointRange& points,
                 const PolygonRange& polygons,
                 const CGAL_NP_CLASS& np)
{
  typedef typename boost::range_value<PolygonRange>::type                   Poly;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, CGAL_NP_CLASS>::type   PointMap;
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));

  if(!os.good())
    return false;

  set_ascii_mode(os); // GOCAD is ASCII only

  set_stream_precision_from_NP(os, np);

  os << "GOCAD TSurf 1\n"
        "HEADER {\n"
        "name:";
  os << fname << "\n";
  os << "*border:on\n"
        "*border*bstone:on\n"
        "}\n"
        "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n"
        "NAME Default\n"
        "AXIS_NAME \"X\" \"Y\" \"Z\"\n"
        "AXIS_UNIT \"m\" \"m\" \"m\"\n"
        "ZPOSITIVE Elevation\n"
        "END_ORIGINAL_COORDINATE_SYSTEM\n"
        "TFACE\n";

  typedef typename Kernel_traits<typename boost::property_traits<PointMap>::value_type>::type K;
  typedef Simple_cartesian<double> SC;
  Cartesian_converter<K,SC> conv;
  for(std::size_t i=0, end=points.size(); i<end; ++i)
    os << "VRTX " << i << " " << conv(get(point_map, points[i])) << "\n";

  for(const Poly& poly : polygons)
  {
    os << "TRGL";
    for(const auto& id : poly)
      os << " " << id;
    os<< "\n";
  }

  os << "END" << std::endl;

  return os.good();
}

} // namespace internal

/*!
 * \ingroup PkgStreamSupportIoFuncsGOCAD
 *
 * \brief writes the content of `points` and `polygons` in `os`, using the \ref IOStreamGocad.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 * \tparam PolygonRange a model of the concept `SequenceContainer` whose `value_type` is itself
 *                      a model of the concept `SequenceContainer` whose `value_type` is
 *                      an unsigned integer type convertible to `std::size_t`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param os the output stream
 * \param points points of the soup of polygons
 * \param polygons a range of polygons. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of the range `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
 *                    of the iterator of `PointRange` and value type is a model of the concept `Kernel::Point_3`}
 *     \cgalParamDefault{`CGAL::Identity_property_map<std::iterator_traits<PointRange::iterator>::value_type>`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e., how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{the precision of the stream `os`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if the writing was successful, `false` otherwise.
*/
template <typename PointRange,
          typename PolygonRange,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(std::ostream& os,
                 const PointRange& points,
                 const PolygonRange& polygons,
                 const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                 , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
                 )
{
  return internal::write_GOCAD(os, "anonymous", points, polygons, np);
}

/*!
 * \ingroup PkgStreamSupportIoFuncsGOCAD
 *
 * \brief writes the content of `points` and `polygons` in `fname`, using the \ref IOStreamGocad.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 * \tparam PolygonRange a model of the concept `SequenceContainer` whose `value_type` is itself
 *                      a model of the concept `SequenceContainer` whose `value_type` is
 *                      an unsigned integer type convertible to `std::size_t`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the path to the output file
 * \param points points of the soup of polygons
 * \param polygons a range of polygons. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of the range `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
 *                    of the iterator of `PointRange` and value type is a model of the concept `Point_3`}
 *     \cgalParamDefault{`CGAL::Identity_property_map<std::iterator_traits<PointRange::iterator>::value_type>`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{`6`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if the writing was successful, `false` otherwise.
*/
template <typename PointRange,
          typename PolygonRange,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(const std::string& fname,
                 const PointRange& points,
                 const PolygonRange& polygons,
                 const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                 , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
                 )
{
  std::ofstream os(fname);
  CGAL::IO::set_mode(os, CGAL::IO::ASCII);
  return internal::write_GOCAD(os, fname.c_str(), points, polygons, np);
}

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_GOCAD_H
