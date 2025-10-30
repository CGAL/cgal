// Copyright (c) 2015-2020  Geometry Factory
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

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/property_map.h>
#include <CGAL/iterator.h>

#include <boost/range/value_type.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <type_traits>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {
namespace internal {

// HEdgesRange = range of std::pair<unsigned int, unsigned int>
// HUVRange = range of std::pair<float, float>
template <class PointRange, class PolygonRange, class ColorOutputIterator, class HEdgesOutputIterator, class HUVOutputIterator>
bool read_PLY(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              HEdgesOutputIterator hedges_out,
              ColorOutputIterator fc_out,
              ColorOutputIterator vc_out,
              HUVOutputIterator huvs_out,
              const bool verbose = false,
              std::enable_if_t<CGAL::is_iterator<ColorOutputIterator>::value>* = nullptr)
{
  typedef typename boost::range_value<PointRange>::type     Point_3;
  typedef CGAL::IO::Color                                   Color_rgb;

  if(!is.good())
  {
    if(verbose)
      std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  internal::PLY_reader reader(verbose);

  if(!(reader.init(is)))
  {
    is.setstate(std::ios::failbit);
    return false;
  }

  for(std::size_t i=0; i<reader.number_of_elements(); ++i)
  {
    internal::PLY_element& element = reader.element(i);

    if(element.name() == "vertex" || element.name() == "vertices")
    {
      bool has_colors = false;
      std::string rtag = "r", gtag = "g", btag = "b";

      if((element.has_property<std::uint8_t>("red") || element.has_property<std::uint8_t>("r")) &&
         (element.has_property<std::uint8_t>("green") || element.has_property<std::uint8_t>("g")) &&
         (element.has_property<std::uint8_t>("blue") || element.has_property<std::uint8_t>("b")))
      {
        has_colors = true;

        if(element.has_property<std::uint8_t>("red"))
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
          internal::PLY_read_number* property = element.property(k);
          property->get(is);

          if(is.fail())
            return false;
        }

        std::tuple<Point_3, std::uint8_t, std::uint8_t, std::uint8_t> new_vertex;
        if(has_colors)
        {
          internal::process_properties(element, new_vertex,
                                       make_ply_point_reader(CGAL::make_nth_of_tuple_property_map<0>(new_vertex)),
                                       std::make_pair(CGAL::make_nth_of_tuple_property_map<1>(new_vertex),
                                                      PLY_property<std::uint8_t>(rtag.c_str())),
                                       std::make_pair(CGAL::make_nth_of_tuple_property_map<2>(new_vertex),
                                                      PLY_property<std::uint8_t>(gtag.c_str())),
                                       std::make_pair(CGAL::make_nth_of_tuple_property_map<3>(new_vertex),
                                                      PLY_property<std::uint8_t>(btag.c_str())));

          *vc_out++ = Color_rgb(get<1>(new_vertex), get<2>(new_vertex), get<3>(new_vertex));
        }
        else
        {
          internal::process_properties(element, new_vertex,
                                       make_ply_point_reader(CGAL::make_nth_of_tuple_property_map<0>(new_vertex)));
        }

        points.push_back(get<0>(new_vertex));
      }
    }
    else if(element.name() == "face" || element.name() == "faces")
    {
      if(element.has_property<std::vector<std::int32_t> >("vertex_indices"))
      {
        internal::read_PLY_faces<std::int32_t>(is, element, polygons, fc_out, "vertex_indices");
      }
      else if(element.has_property<std::vector<std::uint32_t> >("vertex_indices"))
      {
        internal::read_PLY_faces<std::uint32_t>(is, element, polygons, fc_out, "vertex_indices");
      }
      else if(element.has_property<std::vector<std::int32_t> >("vertex_index"))
      {
        internal::read_PLY_faces<std::int32_t>(is, element, polygons, fc_out, "vertex_index");
      }
      else if(element.has_property<std::vector<std::uint32_t> >("vertex_index"))
      {
        internal::read_PLY_faces<std::uint32_t>(is, element, polygons, fc_out, "vertex_index");
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
          internal::PLY_read_number* property = element.property(k);
          property->get(is);

          if(is.fail())
            return false;
        }

        if(has_uv)
        {
          internal::process_properties(element, new_hedge,
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
          internal::process_properties(element, new_hedge,
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
          internal::PLY_read_number* property = element.property(k);
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

/// \cond SKIP_IN_MANUAL

template <class PointRange, class PolygonRange, class ColorRange, class HEdgesRange, class HUVRange>
bool read_PLY(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              HEdgesRange& hedges,
              ColorRange& fcolors,
              ColorRange& vcolors,
              HUVRange& huvs,
              const bool verbose = false,
              std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr)
{
  return internal::read_PLY(is, points, polygons,
                            std::back_inserter(hedges),
                            std::back_inserter(fcolors),
                            std::back_inserter(vcolors),
                            std::back_inserter(huvs),
                            verbose);
}

template <class PointRange, class PolygonRange, class ColorRange>
bool read_PLY(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              ColorRange& fcolors,
              ColorRange& vcolors,
              const bool verbose = false)
{
  std::vector<std::pair<unsigned int, unsigned int> > dummy_pui;
  std::vector<std::pair<float, float> > dummy_pf;

  return internal::read_PLY(is, points, polygons,
                            std::back_inserter(dummy_pui),
                            std::back_inserter(fcolors),
                            std::back_inserter(vcolors),
                            std::back_inserter(dummy_pf),
                            verbose);
}

/// \endcond

/*!
 * \ingroup PkgStreamSupportIoFuncsPLY
 *
 * \brief reads the content of `is` into `points` and `polygons`, using the \ref IOStreamPLY.
 *
 * \attention The polygon soup is not cleared, and the data from the stream are appended.
 *
 * \attention To read a binary file, the flag `std::ios::binary` must be set during the creation of the `ifstream`.
 *
 * \tparam PointRange a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                    whose value type is the point type
 * \tparam PolygonRange a model of the concepts `SequenceContainer` and `BackInsertionSequence`
 *                      whose `value_type` is itself a model of the concept `SequenceContainer`
 *                      and `BackInsertionSequence` whose `value_type` is an unsigned integer type
 *                      convertible to `std::size_t`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param is the input stream
 * \param points points of the soup of polygons
 * \param polygons a range of polygons. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{use_binary_mode}
 *     \cgalParamDescription{indicates whether data should be read in binary (`true`) or in \ascii (`false`)}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{indicates whether output warnings and error messages should be printed or not.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the reading was successful, `false` otherwise.
 */
template <class PointRange, class PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
              , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
              )
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  std::vector<std::pair<unsigned int, unsigned int> > dummy_pui;
  std::vector<std::pair<float, float> > dummy_pf;

  return internal::read_PLY(is, points, polygons, std::back_inserter(dummy_pui),
                            choose_parameter(get_parameter(np, internal_np::face_color_output_iterator),
                                             CGAL::Emptyset_iterator()),
                            choose_parameter(get_parameter(np, internal_np::vertex_color_output_iterator),
                                             CGAL::Emptyset_iterator()),
                            std::back_inserter(dummy_pf),
                            choose_parameter(get_parameter(np, internal_np::verbose), true));
}

/*!
 * \ingroup PkgStreamSupportIoFuncsPLY
 *
 * \brief reads the content of `fname` into `points` and `polygons`, using the \ref IOStreamPLY.
 *
 * \attention The polygon soup is not cleared, and the data from the file are appended.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * \tparam PolygonRange a model of the concepts `SequenceContainer` and `BackInsertionSequence`
 *                      whose `value_type` is itself a model of the concept `SequenceContainer`
 *                      and `BackInsertionSequence` whose `value_type` is an integer type
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the path to the input file
 * \param points points of the soup of polygons
 * \param polygons a range of polygons. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{use_binary_mode}
 *     \cgalParamDescription{indicates whether data should be read in binary (`true`) or in \ascii (`false`)}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{indicates whether output warnings and error messages should be printed or not.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the reading was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const std::string& fname,
              PointRange& points,
              PolygonRange& polygons,
              const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
              , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
              )
{
  const bool binary = parameters::choose_parameter(parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ifstream is(fname, std::ios::binary);
    CGAL::IO::set_mode(is, CGAL::IO::BINARY);
    return read_PLY(is, points, polygons, np);
  }
  else
  {
    std::ifstream is(fname);
    CGAL::IO::set_mode(is, CGAL::IO::ASCII);
    return read_PLY(is, points, polygons, np);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

// @todo comments

/*!
 * \ingroup PkgStreamSupportIoFuncsPLY
 *
 * \brief writes the content of `points` and `polygons` in `out`, using the \ref IOStreamPLY.
 *
 * \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation
 *            of the `ofstream`, and the \link PkgStreamSupportEnumRef `IO::Mode` \endlink
 *            of the stream must be set to `BINARY`.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 * \tparam PolygonRange a model of the concept `SequenceContainer`
 *                      whose `value_type` is itself a model of the concept `SequenceContainer`
 *                      whose `value_type` is an unsigned integer type convertible to `std::size_t`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param out the output stream
 * \param points points of the soup of polygons
 * \param polygons a range of polygons. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{the precision of the stream `os`}
 *     \cgalParamExtra{This parameter is only meaningful while using \ascii encoding.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if the writing was successful, `false` otherwise.
 */
template <class PointRange, class PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS >
bool write_PLY(std::ostream& out,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
               )
{
  typedef typename boost::range_value<PointRange>::type Point_3;
  typedef typename boost::range_value<PolygonRange>::type Polygon_3;

  if(!out.good())
    return false;

  set_stream_precision_from_NP(out, np);

  // Write header
  out << "ply" << std::endl
      << ((get_mode(out) == BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
      << "comment Generated by the CGAL library" << std::endl
      << "element vertex " << points.size() << std::endl;

  internal::output_property_header(out, make_ply_point_writer(CGAL::Identity_property_map<Point_3>()));

  out << "element face " << polygons.size() << std::endl;

  internal::output_property_header(out, std::make_pair(CGAL::Identity_property_map<Polygon_3>(),
                                                       PLY_property<std::vector<int> >("vertex_indices")));

  out << "end_header" << std::endl;

  for(std::size_t i=0; i<points.size(); ++i)
    internal::output_properties(out, points.begin() + i,
                                make_ply_point_writer(CGAL::Identity_property_map<Point_3>()));

  for(std::size_t i=0; i<polygons.size(); ++i)
    internal::output_properties(out, polygons.begin() + i,
                                std::make_pair(CGAL::Identity_property_map<Polygon_3>(),
                                               PLY_property<std::vector<int> >("vertex_indices")));

  return out.good();
}

/*!
 * \ingroup PkgStreamSupportIoFuncsPLY
 *
 * \brief writes the content of `points` and `polygons` in  the file `fname`, using the \ref IOStreamPLY.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 * \tparam PolygonRange a model of the concept `SequenceContainer`
 *                      whose `value_type` is itself a model of the concept `SequenceContainer`
 *                      whose `value_type` is an unsigned integer type convertible to `std::size_t`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the path to the output file
 * \param points points of the soup of polygons
 * \param polygons a range of polygons. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{use_binary_mode}
 *     \cgalParamDescription{indicates whether data should be written in binary (`true`) or in \ascii (`false`)}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{`6`}
 *     \cgalParamExtra{This parameter is only meaningful while using \ascii encoding.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if the writing was successful, `false` otherwise.
 */
template <class PointRange, class PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS >
bool write_PLY(const std::string& fname,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif  // DOXYGEN_RUNNING
               )
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ofstream os(fname, std::ios::binary);
    CGAL::IO::set_mode(os, CGAL::IO::BINARY);
    return write_PLY(os, points, polygons, np);
  }
  else
  {
    std::ofstream os(fname);
    CGAL::IO::set_mode(os, CGAL::IO::ASCII);
    return write_PLY(os, points, polygons, np);
  }
}



#ifdef DOXYGEN_RUNNING
/**
   \ingroup PkgStreamSupportIoFuncsPLY

   Class used to identify a %PLY property as a type and a name.

   \sa `read_PLY_with_properties()`
*/
template <typename T>
struct PLY_property
{
  typedef T type;
  const char* name;
  PLY_property(const char* name) : name(name) { }
};

/**
   \ingroup PkgStreamSupportIoFuncsPLY

   Generates a %PLY property handler to read 3D points. Points are
   constructed from the input using 3 %PLY properties of type `FT`
   and named `x`, `y` and `z`. `FT` is `float` if the points use
   `CGAL::Simple_cartesian<float>` and `double` otherwise.

   \tparam PointMap the property map used to store points.

   \sa `read_PLY_with_properties()`
   \sa \ref IOStreamPLY
*/
template <typename PointMap>
std::tuple<PointMap,
           typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3,
           PLY_property<FT>, PLY_property<FT>, PLY_property<FT> >
make_ply_point_reader(PointMap point_map);

/**
   \ingroup PkgStreamSupportIoFuncsPLY

   Generates a %PLY property handler to read 3D normal
   vectors. Vectors are constructed from the input using 3 PLY
   properties of type `FT` and named `nx`, `ny` and `nz`. `FT`
   is `float` if the points use `CGAL::Simple_cartesian<float>` and
   `double` otherwise.

   \tparam VectorMap the property map used to store vectors.

   \sa `read_PLY_with_properties()`
   \sa \ref IOStreamPLY
*/
template <typename VectorMap>
std::tuple<VectorMap,
           typename Kernel_traits<typename VectorMap::value_type>::Kernel::Construct_vector_3,
           PLY_property<FT>, PLY_property<FT>, PLY_property<FT> >
make_ply_normal_reader(VectorMap normal_map);

#endif // DOXYGEN_RUNNING

/**
  \ingroup PkgStreamSupportIoFuncsPLY

  \brief reads user-selected points properties from a .ply stream (ASCII or binary).

  Potential additional point properties and faces are ignored.

  Properties are handled through a variadic list of property
  handlers. A `PropertyHandler` can either be:

  - A `std::pair<PropertyMap, PLY_property<T> >` if the user wants
  to read a %PLY property as a scalar value T (for example, storing
  an `int` %PLY property into an `int` variable).

  - A `std::tuple<PropertyMap, Constructor,
  PLY_property<T>...>` if the user wants to use one or several PLY
  properties to construct a complex object (for example, storing 3
  `uchar` %PLY properties into a %Color object that can for example
  be a `std::array<unsigned char, 3>`). In that case, the
  second element of the tuple should be a functor that constructs
  the value type of `PropertyMap` from N objects of types `T`.

  \attention To read a binary file, the flag `std::ios::binary` must be set during the creation of the `ifstream`.

  \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
  It must be a model of `DefaultConstructible` and defaults to `value_type_traits<PointOutputIterator>::%type`.
  It can be omitted when the default is fine.
  \tparam PointOutputIterator iterator over output points.
  \tparam PropertyHandler handlers to recover properties.

  \returns `true` if reading was successful, `false` otherwise.

  \sa \ref IOStreamPLY
  \sa `make_ply_point_reader()`
  \sa `make_ply_normal_reader()`
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename ... PropertyHandler>
bool read_PLY_with_properties(std::istream& is,
                              PointOutputIterator output,
                              PropertyHandler&& ... properties);



/**
   \ingroup PkgStreamSupportIoFuncsPLY

   \brief reads points (positions + normals, if available), using the \ref IOStreamPLY.

   Potential additional point properties and faces are ignored.

   \tparam OutputIteratorValueType type of objects that can be put in `PointOutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<PointOutputIterator>::%type`.
   It can be omitted when the default is fine.
   \tparam PointOutputIterator iterator over output points.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param fname input file name.
   \param output output iterator over points.
   \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamNBegin{use_binary_mode}
       \cgalParamDescription{indicates whether data should be read in binary (`true`) or in \ascii (`false`)}
       \cgalParamType{Boolean}
       \cgalParamDefault{`true`}
     \cgalParamNEnd

     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals in the input stream are ignored.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if reading was successful, `false` otherwise.

   \sa \ref IOStreamPLY
   \sa `read_PLY_with_properties()`
*/
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const std::string& fname,
              PointOutputIterator output,
              const CGAL_NP_CLASS& np = parameters::default_values()
              #ifndef DOXYGEN_RUNNING
              , std::enable_if_t<CGAL::is_iterator<PointOutputIterator>::value>* = nullptr
              #endif
              );



#ifdef DOXYGEN_RUNNING // Document some parts from Stream_support here for convenience
  /**
     \ingroup PkgStreamSupportIoFuncsPLY

     Generates a %PLY property handler to write 3D points. Points are
     written as 3 %PLY properties of type `FT` and named `x`, `y` and
     `z`. `FT` is `float` if the points use
     `CGAL::Simple_cartesian<float>` and `double` otherwise.

     \tparam PointMap the property map used to store points.

     \sa `write_PLY_with_properties()`
     \sa \ref IOStreamPLY
  */
  template <typename PointMap>
  std::tuple<PointMap, PLY_property<FT>, PLY_property<FT>, PLY_property<FT> >
  make_ply_point_writer(PointMap point_map);

  /**
     \ingroup PkgStreamSupportIoFuncsPLY

     Generates a %PLY property handler to write 3D normal
     vectors. Vectors are written as 3 %PLY properties of type `FT`
     and named `nx`, `ny` and `nz`. `FT` is `float` if the vectors use
     `CGAL::Simple_cartesian<float>` and `double` otherwise.

     \tparam VectorMap the property map used to store vectors.

     \sa `write_PLY_with_properties()`
     \sa \ref IOStreamPLY
  */
  template <typename VectorMap>
  std::tuple<VectorMap, PLY_property<FT>, PLY_property<FT>, PLY_property<FT> >
  make_ply_normal_writer(VectorMap normal_map);
#endif


/**
   \ingroup PkgStreamSupportIoFuncsPLY

   \brief writes the range of `points` with properties using \ref IOStreamPLY.

   Properties are handled through a variadic list of property
   handlers. A `PropertyHandler` can either be:

   - A `std::pair<PropertyMap, PLY_property<T> >` if the user wants
   to write a scalar value T as a %PLY property (for example, writing
   an `int` variable as an `int` %PLY property).

   - A `std::tuple<PropertyMap, PLY_property<T>...>` if the
   user wants to write a complex object as several %PLY
   properties. In that case, a specialization of `Output_rep` must
   be provided for `PropertyMap::value_type` that handles both ASCII
   and binary output (see `CGAL::IO::get_mode()`).

   \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation
              of the `ofstream`, and the \link PkgStreamSupportEnumRef `IO::Mode` \endlink
              of the stream must be set to `BINARY`.

   \tparam PointRange is a model of `ConstRange`. The value type of
                      its iterator is the key type of the `PropertyMap` objects provided
                      within the `PropertyHandler` parameter.
   \tparam PropertyHandler handlers to recover properties.

   \returns `true` if writing was successful, `false` otherwise.

   \sa \ref IOStreamPLY
   \sa `make_ply_point_writer()`
   \sa `make_ply_normal_writer()`
*/
template <typename PointRange,
          typename ... PropertyHandler>
  bool write_PLY_with_properties(std::ostream& os, ///< output stream.
                                 const PointRange& points, ///< input point range.
                                 PropertyHandler&& ... properties); ///< parameter pack of property handlers



/**
   \ingroup PkgStreamSupportIoFuncsPLY

   \brief writes the range of `points` (positions + normals, if available) using \ref IOStreamPLY.

   \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation
              of the `ofstream`, and the \link PkgStreamSupportEnumRef `IO::Mode` \endlink
              of the stream must be set to `BINARY`.

   \tparam PointRange is a model of `ConstRange`. The value type of
                      its iterator is the key type of the named parameter `point_map`.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param os output stream
   \param points input point range
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals are not written in the output stream.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd

     \cgalParamNBegin{stream_precision}
       \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
       \cgalParamType{int}
       \cgalParamDefault{the precision of the stream `os`}
       \cgalParamExtra{This parameter is only meaningful while using \ascii encoding.}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if writing was successful, `false` otherwise.

   \sa `write_PLY_with_properties()`
*/
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(std::ostream& os,
               const PointRange& points,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PointRange>::value>* = nullptr
#endif
               );


/**
   \ingroup PkgStreamSupportIoFuncsPLY

   \brief writes the range of `points` (positions + normals, if available) using \ref IOStreamPLY.

   \tparam PointRange is a model of `ConstRange`. The value type of
                      its iterator is the key type of the named parameter `point_map`.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param filename the path to the output file
   \param points input point range
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{use_binary_mode}
       \cgalParamDescription{indicates whether data should be written in binary (`true`) or in \ascii (`false`)}
       \cgalParamType{Boolean}
       \cgalParamDefault{`true`}
     \cgalParamNEnd

     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals are not written in the output file.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd

     \cgalParamNBegin{stream_precision}
       \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
       \cgalParamType{int}
       \cgalParamDefault{`6`}
       \cgalParamExtra{This parameter is only meaningful while using \ascii encoding.}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if writing was successful, `false` otherwise.

   \sa `write_PLY_with_properties()`
*/
template <typename PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(const std::string& filename,
               const PointRange& points,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PointRange>::value>* = nullptr
#endif
               );

} // namespace IO

} // namespace CGAL

#include <CGAL/IO/PLY/read_ply_points.h>
#include <CGAL/IO/PLY/write_ply_points.h>


#endif // CGAL_IO_PLY_H
