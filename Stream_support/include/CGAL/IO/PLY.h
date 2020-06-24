// Copyright (c) 2015  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_IO_PLY_H
#define CGAL_IO_PLY_H

#include <CGAL/IO/PLY/PLY_reader.h>
#include <CGAL/IO/PLY/PLY_writer.h>
#include <CGAL/IO/helpers.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/property_map.h>
#include <CGAL/iterator.h>

#include <boost/range/value_type.hpp>
#include <boost/utility/enable_if.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {
namespace internal {

// HEdgesRange" = range of std::pair<unsigned int, unsigned int>
// HUVRange = range of std::pair<float, float>
template <class PointRange, class PolygonRange, class ColorOutputIterator, class HEdgesOutputIterator, class HUVOutputIterator>
bool read_PLY(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              HEdgesOutputIterator hedges_out,
              ColorOutputIterator fc_out,
              ColorOutputIterator vc_out,
              HUVOutputIterator huvs_out,
              const bool verbose = true,
              typename std::enable_if<CGAL::is_iterator<ColorOutputIterator>::value>::type* = nullptr)
{
  typedef typename boost::range_value<PointRange>::type     Point_3;
  typedef CGAL::Color                                       Color_rgb;

  if(!is.good())
  {
    if(verbose)
      std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  IO::internal::PLY_reader reader(verbose);

  if(!(reader.init(is)))
  {
    is.setstate(std::ios::failbit);
    return false;
  }

  for(std::size_t i=0; i<reader.number_of_elements(); ++i)
  {
    IO::internal::PLY_element& element = reader.element(i);

    if(element.name() == "vertex" || element.name() == "vertices")
    {
      bool has_colors = false;
      std::string rtag = "r", gtag = "g", btag = "b";

      if((element.has_property<boost::uint8_t>("red") || element.has_property<boost::uint8_t>("r")) &&
         (element.has_property<boost::uint8_t>("green") || element.has_property<boost::uint8_t>("g")) &&
         (element.has_property<boost::uint8_t>("blue") || element.has_property<boost::uint8_t>("b")))
      {
        has_colors = true;

        if(element.has_property<boost::uint8_t>("red"))
        {
          rtag = "red";
          gtag = "green";
          btag = "blue";
        }
      }

      for(std::size_t j=0; j<element.number_of_items(); ++j)
      {
        for(std::size_t k=0; k<element.number_of_properties(); ++k)
        {
          IO::internal::PLY_read_number* property = element.property(k);
          property->get(is);

          if(is.fail())
            return false;
        }

        std::tuple<Point_3, boost::uint8_t, boost::uint8_t, boost::uint8_t> new_vertex;
        if(has_colors)
        {
          IO::internal::process_properties(element, new_vertex,
                                           make_ply_point_reader(CGAL::make_nth_of_tuple_property_map<0>(new_vertex)),
                                           std::make_pair(CGAL::make_nth_of_tuple_property_map<1>(new_vertex),
                                                          PLY_property<boost::uint8_t>(rtag.c_str())),
                                           std::make_pair(CGAL::make_nth_of_tuple_property_map<2>(new_vertex),
                                                          PLY_property<boost::uint8_t>(gtag.c_str())),
                                           std::make_pair(CGAL::make_nth_of_tuple_property_map<3>(new_vertex),
                                                          PLY_property<boost::uint8_t>(btag.c_str())));

          *vc_out++ = Color_rgb(get<1>(new_vertex), get<2>(new_vertex), get<3>(new_vertex));
        }
        else
        {
          IO::internal::process_properties(element, new_vertex,
                                            make_ply_point_reader(CGAL::make_nth_of_tuple_property_map<0>(new_vertex)));
        }

        points.push_back(get<0>(new_vertex));
      }
    }
    else if(element.name() == "face" || element.name() == "faces")
    {
      if(element.has_property<std::vector<boost::int32_t> >("vertex_indices"))
      {
        IO::internal::read_PLY_faces<boost::int32_t>(is, element, polygons, fc_out, "vertex_indices");
      }
      else if(element.has_property<std::vector<boost::uint32_t> >("vertex_indices"))
      {
        IO::internal::read_PLY_faces<boost::uint32_t>(is, element, polygons, fc_out, "vertex_indices");
      }
      else if(element.has_property<std::vector<boost::int32_t> >("vertex_index"))
      {
        IO::internal::read_PLY_faces<boost::int32_t>(is, element, polygons, fc_out, "vertex_index");
      }
      else if(element.has_property<std::vector<boost::uint32_t> >("vertex_index"))
      {
        IO::internal::read_PLY_faces<boost::uint32_t>(is, element, polygons, fc_out, "vertex_index");
      }
      else
      {
        if(verbose)
          std::cerr << "Error: can't find vertex indices in PLY input" << std::endl;
        return false;
      }
    }
    else if(element.name() == "halfedge" )
    {
      bool has_uv = false;
      std::string stag = "source", ttag = "target", utag = "u", vtag = "v";
      if(element.has_property<unsigned int>("source") &&
         element.has_property<unsigned int>("target") &&
         element.has_property<float>("u") &&
         element.has_property<float>("v"))
      {
        has_uv = true;
      }

      std::tuple<unsigned int, unsigned int, float, float, float> new_hedge;
      for(std::size_t j=0; j<element.number_of_items(); ++j)
      {
        for(std::size_t k=0; k<element.number_of_properties(); ++k)
        {
          IO::internal::PLY_read_number* property = element.property(k);
          property->get(is);

          if(is.fail())
            return false;
        }

        if(has_uv)
        {
          IO::internal::process_properties(element, new_hedge,
                                           std::make_pair(CGAL::make_nth_of_tuple_property_map<0>(new_hedge),
                                                          PLY_property<unsigned int>(stag.c_str())),
                                           std::make_pair(CGAL::make_nth_of_tuple_property_map<1>(new_hedge),
                                                          PLY_property<unsigned int>(ttag.c_str())),
                                           std::make_pair(CGAL::make_nth_of_tuple_property_map<2>(new_hedge),
                                                          PLY_property<float>(utag.c_str())),
                                           std::make_pair(CGAL::make_nth_of_tuple_property_map<3>(new_hedge),
                                                          PLY_property<float>(vtag.c_str())));

          *hedges_out++ = std::make_pair(get<0>(new_hedge), get<1>(new_hedge));
          *huvs_out++ = std::make_pair(get<2>(new_hedge), get<3>(new_hedge));
        }
        else
        {
          IO::internal::process_properties(element, new_hedge,
                                           std::make_pair(CGAL::make_nth_of_tuple_property_map<0>(new_hedge),
                                                          PLY_property<unsigned int>(stag.c_str())),
                                           std::make_pair(CGAL::make_nth_of_tuple_property_map<1>(new_hedge),
                                                          PLY_property<unsigned int>(ttag.c_str())));
        }
      }
    }
    else // Read other elements and ignore
    {
      for(std::size_t j=0; j<element.number_of_items(); ++j)
      {
        for(std::size_t k=0; k<element.number_of_properties(); ++k)
        {
          IO::internal::PLY_read_number* property = element.property(k);
          property->get(is);
          if(is.fail())
            return false;
        }
      }
    }
  }
  return !is.fail();
}

} // namespace internal
} // namespace IO

template <class PointRange, class PolygonRange, class ColorRange, class HEdgesRange, class HUVRange>
bool read_PLY(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              HEdgesRange& hedges,
              ColorRange& fcolors,
              ColorRange& vcolors,
              HUVRange& huvs,
              const bool verbose = true,
              typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return IO::internal::read_PLY(is, points, polygons, std::back_inserter(hedges), std::back_inserter(fcolors), std::back_inserter(vcolors), std::back_inserter(huvs), verbose);
}

template <class PointRange, class PolygonRange, class ColorRange>
bool read_PLY(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              ColorRange& fcolors,
              ColorRange& vcolors,
              const bool verbose = true)
{
  std::vector<std::pair<unsigned int, unsigned int> > dummy_pui;
  std::vector<std::pair<float, float> > dummy_pf;

  return IO::internal::read_PLY(is, points, polygons, std::back_inserter(dummy_pui), std::back_inserter(fcolors), std::back_inserter(vcolors), std::back_inserter(dummy_pf), verbose);
}

/*!
 * \ingroup PkgStreamSupportIoFuncsPLY
 *
 * \brief Reads the content of `is` into `points` and `polygons`, using the \ref IOStreamPLY.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * \tparam PolygonRange a model of the concept `SequenceContainer`
 *                      whose value_type is itself a model of the concept `SequenceContainer`
 *                      whose value_type is an integer type.
 *
 * \param is the input stream
 * \param points points of the soup of polygons.
 * \param polygons a `PolygonRange`. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param verbose: if `true`, will output warnings and error messages.
 *
 * \attention Be mindful of the flag `std::ios::binary` flag when creating the `ifstream` when reading a binary file
 *
 * \returns `true` if the reading was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              const CGAL_BGL_NP_CLASS& np,
              const bool verbose = true
#ifndef DOXYGEN_RUNNING
              , typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr
#endif
              )
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  std::vector<std::pair<unsigned int, unsigned int> > dummy_pui;
  std::vector<std::pair<float, float> > dummy_pf;

  return IO::internal::read_PLY(is, points, polygons, std::back_inserter(dummy_pui),
                                choose_parameter(get_parameter(np, internal_np::face_color_output_iterator),
                                                 CGAL::Emptyset_iterator()),
                                choose_parameter(get_parameter(np, internal_np::vertex_color_output_iterator),
                                                 CGAL::Emptyset_iterator()),
                                std::back_inserter(dummy_pf), verbose);
}


template <class PointRange, class PolygonRange>
bool read_PLY(std::istream& is, PointRange& points, PolygonRange& polygons,
              const bool verbose = true, typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return read_PLY(is, points, polygons, parameters::all_default(), verbose);
}

/*!
 * \ingroup PkgStreamSupportIoFuncsPLY
 *
 * \brief Reads the content of `fname` into `points` and `polygons`, using the \ref IOStreamPLY.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * \tparam PolygonRange a model of the concept `SequenceContainer`
 *                      whose value_type is itself a model of the concept `SequenceContainer`
 *                      whose value_type is an integer type
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the path to the input file
 * \param points points of the soup of polygons.
 * \param polygons a `PolygonRange`. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{use_binary_mode}
 *     \cgalParamDescription{indicates whether data should be read in binary (`true`) or in ASCII (`false`)}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the reading was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const char* fname,
              PointRange& points,
              PolygonRange& polygons,
              const CGAL_BGL_NP_CLASS& np,
              const bool verbose = true,
#ifndef DOXYGEN_RUNNING
              typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr
#endif
              )
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ifstream is(fname, std::ios::binary);
    CGAL::set_mode(is, CGAL::IO::BINARY);
    return read_PLY(is, points, polygons, np, verbose);
  }
  else
  {
    std::ifstream is(fname);
    CGAL::set_mode(is, CGAL::IO::ASCII);
    return read_PLY(is, points, polygons, np, verbose);
  }
}

template <typename PointRange, typename PolygonRange>
bool read_PLY(const char* fname, PointRange& points, PolygonRange& polygons,
              const bool verbose = true, typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return read_PLY(fname, points, polygons, parameters::all_default(), verbose);
}

template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const std::string& fname, PointRange& points, PolygonRange& polygons, const CGAL_BGL_NP_CLASS& np,
              const bool verbose = true, typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return read_PLY(fname.c_str(), points, polygons, np, verbose);
}

template <typename PointRange, typename PolygonRange>
bool read_PLY(const std::string& fname, PointRange& points, PolygonRange& polygons,
              const bool verbose = true, typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return read_PLY(fname.c_str(), points, polygons, parameters::all_default(), verbose);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

// @todo comments

/*!
 * \ingroup PkgStreamSupportIoFuncsPLY
 *
 * writes the content of `points` and `polygons` in `out`, using the \ref IOStreamPLY.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * \tparam PolygonRange a model of the concept `SequenceContainer`
 *                      whose value_type is itself a model of the concept `SequenceContainer`
 *                      whose value_type is an integer type.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param out the output stream
 * \param points points of the soup of polygons.
 * \param polygons a `PolygonRange`. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{`6`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \attention Be mindful of the flag `std::ios::binary` flag when creating the `ofstream` when writing a binary file
 *
 * \return `true` if the writing was successful, `false` otherwise.
 */
template <class PointRange, class PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS >
bool write_PLY(std::ostream& out,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr
#endif
               )
{
  typedef typename boost::range_value<PointRange>::type Point_3;
  typedef typename boost::range_value<PolygonRange>::type Polygon_3;

  if(!out.good())
    return false;

  const int precision = parameters::choose_parameter(parameters::get_parameter(np, internal_np::stream_precision), 6);
  out << std::setprecision(precision);

  // Write header
  out << "ply" << std::endl
      << ((get_mode(out) == IO::BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
      << "comment Generated by the CGAL library" << std::endl
      << "element vertex " << points.size() << std::endl;

  IO::internal::output_property_header(out, make_ply_point_writer(CGAL::Identity_property_map<Point_3>()));

  out << "element face " << polygons.size() << std::endl;

  IO::internal::output_property_header(out, std::make_pair(CGAL::Identity_property_map<Polygon_3>(),
                                                            PLY_property<std::vector<int> >("vertex_indices")));

  out << "end_header" << std::endl;

  for(std::size_t i=0; i<points.size(); ++i)
    IO::internal::output_properties(out, points.begin() + i,
                                    make_ply_point_writer(CGAL::Identity_property_map<Point_3>()));

  for(std::size_t i=0; i<polygons.size(); ++i)
    IO::internal::output_properties(out, polygons.begin() + i,
                                    std::make_pair(CGAL::Identity_property_map<Polygon_3>(),
                                                   PLY_property<std::vector<int> >("vertex_indices")));

  return out.good();
}

template <class PointRange, class PolygonRange>
bool write_PLY(std::ostream& out, const PointRange& points, const PolygonRange& polygons,
               typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return write_PLY(out, points, polygons, parameters::all_default());
}

/*!
 * \ingroup PkgStreamSupportIoFuncsPLY
 *
 * \brief Writes the content of `points` and `polygons` in  the file `fname`, using the \ref IOStreamPLY.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * \tparam PolygonRange a model of the concept `SequenceContainer`
 *                      whose value_type is itself a model of the concept `SequenceContainer`
 *                      whose value_type is an integer type.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the path to the output file
 * \param points points of the soup of polygons.
 * \param polygons a `PolygonRange`. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{use_binary_mode}
 *     \cgalParamDescription{indicates whether data should be written in binary (`true`) or in ASCII (`false`)}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{`6`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if the writing was successful, `false` otherwise.
 */
template <class PointRange, class PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS >
bool write_PLY(const char* fname,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr
#endif
               )
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ofstream os(fname, std::ios::binary);
    CGAL::set_mode(os, CGAL::IO::BINARY);
    return write_PLY(os, points, polygons, np);
  }
  else
  {
    std::ofstream os(fname);
    CGAL::set_mode(os, CGAL::IO::ASCII);
    return write_PLY(os, points, polygons, np);
  }
}

template <class PointRange, class PolygonRange>
bool write_PLY(const char* fname, const PointRange& points, const PolygonRange& polygons,
               typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return write_PLY(fname, points, polygons, parameters::all_default());
}

template <class PointRange, class PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS >
bool write_PLY(const std::string& fname, const PointRange& points, const PolygonRange& polygons, const CGAL_BGL_NP_CLASS& np,
               typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return write_PLY(fname.c_str(), points, polygons, np);
}

template <class PointRange, class PolygonRange>
bool write_PLY(const std::string& fname, const PointRange& points, const PolygonRange& polygons,
               typename boost::enable_if<IO::internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return write_PLY(fname.c_str(), points, polygons, parameters::all_default());
}

} // namespace CGAL

#endif // CGAL_IO_PLY_H
