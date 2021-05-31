// Copyright (c) 2015-2020  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Lutz Kettner
//             Andreas Fabri
//             Maxime Gimeno

#ifndef CGAL_IO_OBJ_H
#define CGAL_IO_OBJ_H

#include <CGAL/IO/OBJ/File_writer_wavefront.h>
#include <CGAL/IO/Generic_writer.h>
#include <CGAL/IO/io.h>
#include <CGAL/IO/helpers.h>

#include <CGAL/Container_helper.h>

#include <boost/range/value_type.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/boost/graph/Named_function_parameters.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#endif

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {
namespace internal {

template <typename PointRange, typename PolygonRange, typename VertexNormalOutputIterator, typename VertexTextureOutputIterator>
bool read_OBJ(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              VertexNormalOutputIterator,
              VertexTextureOutputIterator,
              const bool verbose = false)
{
  if(!is.good())
  {
    if(verbose)
      std::cerr<<"File doesn't exist."<<std::endl;
    return false;
  }

  typedef typename boost::range_value<PointRange>::type                               Point;

  set_ascii_mode(is); // obj is ASCII only

  int mini(1), maxi(-1);
  std::string s;
  Point p;

  std::string line;
  bool tex_found(false), norm_found(false);
  while(getline(is, line))
  {
    if(line.empty())
      continue;

    std::istringstream iss(line);
    if(!(iss >> s))
      continue; // can't read anything on the line, whitespace only?

    if(s == "v")
    {
      if(!(iss >> p))
      {
        if(verbose)
          std::cerr << "error while reading OBJ vertex." << std::endl;
        return false;
      }

      points.push_back(p);
    }
    else if(s == "vt")
    {
      tex_found = true;
    }
    else if(s == "vn")
    {
      norm_found = true;
    }
    else if(s == "f")
    {
      int i;
      polygons.emplace_back();
      while(iss >> i)
      {
        if(i < 1)
        {
          const std::size_t n = polygons.back().size();
          ::CGAL::internal::resize(polygons.back(), n + 1);
          polygons.back()[n] = static_cast<int>(points.size()) + i; // negative indices are relative references
          if(i < mini)
            mini = i;
        }
        else
        {
          const std::size_t n = polygons.back().size();
          ::CGAL::internal::resize(polygons.back(), n + 1);
          polygons.back()[n] = i - 1;
          if(i-1 > maxi)
            maxi = i-1;
        }

        // the format can be "f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 ..." and we only read vertex ids for now,
        // so skip to the next vertex
        iss.ignore(256, ' ');
      }

      if(iss.bad())
        return false;
    }
    else if(s.front() == '#')
    {
      // this is a commented line, ignored
    }
    else if(s == "vp" ||
            // Display
            s == "bevel" || s == "lod" || s == "ctech" || s == "c_interp" || s == "usemap" || s == "usemtl" ||
            s == "stech" || s == "d_interp" || s == "mtllib" || s == "shadow_obj" || s == "trace_obj" ||
            // groups
            s == "o" || s == "g" || s == "s" ||
            // Free
            s == "p" || s == "cstype" || s == "deg" || s == "step" || s == "bmat" || s == "con" ||
            s == "curv" || s == "curv2" || s == "surf" || s == "parm" || s == "trim" || s == "hole" ||
            s == "scrv" || s == "sp" || s == "end" ||
            s == "con" || s == "surf_1" || s == "q0_1" || s == "q1_1" || s == "curv2d_1" ||
            s == "surf_2" || s == "q0_2" || s == "q1_2" || s == "curv2d_2" ||
            // supersed statements
            s == "bsp" || s == "bzp" || s == "cdc" || s == "cdp" || s == "res")
    {
      // valid, but unsupported
    }
    else
    {
      if(verbose)
        std::cerr << "error: unrecognized line: " << s << std::endl;
      return false;
    }
  }

  if(norm_found && verbose)
    std::cout<<"NOTE: normals were found in this file, but were discarded."<<std::endl;
  if(tex_found && verbose)
    std::cout<<"NOTE: textures were found in this file, but were discarded."<<std::endl;

  if(points.empty() || polygons.empty())
  {
    if(verbose)
      std::cerr << "warning: empty file?" << std::endl;
    return false;
  }

  if(maxi > static_cast<int>(points.size()) || mini < -static_cast<int>(points.size()))
  {
    if(verbose)
      std::cerr << "error: invalid face index" << std::endl;
    return false;
  }

  return !is.bad();
}

} // namespace internal

/// \ingroup PkgStreamSupportIoFuncsOBJ
///
/// \brief reads the content of `is` into `points` and `polygons`, using the \ref IOStreamOBJ.
///
/// \attention The polygon soup is not cleared, and the data from the stream are appended.
///
/// \tparam PointRange a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
///                    whose value type is the point type
/// \tparam PolygonRange a model of the concepts `SequenceContainer` and `BackInsertionSequence`
///                      whose `value_type` is itself a model of the concepts `SequenceContainer`
///                      and `BackInsertionSequence` whose `value_type` is an unsigned integer type
///                      convertible to `std::size_t`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param is the input stream
/// \param points points of the soup of polygons
/// \param polygons a range of polygons. Each element in it describes a polygon
///        using the indices of the points in `points`.
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{verbose}
///     \cgalParamDescription{indicates whether output warnings and error messages should be printed or not.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`false`}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \returns `true` if the reading was successful, `false` otherwise.
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
              , typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr
#endif
              )
{
  const bool verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);

  return internal::read_OBJ(is, points, polygons,
                            CGAL::Emptyset_iterator(), CGAL::Emptyset_iterator(),
                            verbose);
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange>
bool read_OBJ(std::istream& is, PointRange& points, PolygonRange& polygons,
              typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return read_OBJ(is, points, polygons, parameters::all_default());
}

/// \endcond

/// \ingroup PkgStreamSupportIoFuncsOBJ
///
/// \brief reads the content of the file `fname` into `points` and `polygons`, using the \ref IOStreamOBJ.
///
/// \attention The polygon soup is not cleared, and the data from the file are appended.
///
/// \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
/// \tparam PolygonRange a model of the concepts `SequenceContainer` and `BackInsertionSequence`
///                      whose `value_type` is itself a model of the concept `SequenceContainer`
///                      and `BackInsertionSequence` whose `value_type` is an unsigned integer type
///                      convertible to `std::size_t`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param fname the path to the input file
/// \param points points of the soup of polygons
/// \param polygons a range of polygons. Each element in it describes a polygon
///        using the indices of the points in `points`.
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{verbose}
///     \cgalParamDescription{indicates whether output warnings and error messages should be printed or not.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`false`}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \returns `true` if the reading was successful, `false` otherwise.
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(const std::string& fname,
              PointRange& points,
              PolygonRange& polygons,
              const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
              , typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr
#endif
              )
{
  std::ifstream is(fname);
  CGAL::IO::set_mode(is, CGAL::IO::ASCII);
  return read_OBJ(is, points, polygons, np);
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange>
bool read_OBJ(const std::string& fname, PointRange& points, PolygonRange& polygons,
              typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return read_OBJ(fname, points, polygons, parameters::all_default());
}

/// \endcond

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 * \ingroup PkgStreamSupportIoFuncsOBJ
 *
 * \brief writes the content of `points` and `polygons` in `os`, using the \ref IOStreamOBJ.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 * \tparam PolygonRange a model of the concept `SequenceContainer`
 *                      whose `value_type` is itself a model of the concept `SequenceContainer`
 *                      whose `value_type` is an unsigned integer type convertible to `std::size_t`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param os the output stream
 * \param points points of the soup of polygons
 * \param polygons a range of polygons. Each element in it describes a polygon
 *        using the indices of the points in `points`.
 * \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{`the precision of the stream `os``}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if the writing was successful, `false` otherwise.
 */
template <typename PointRange,
          typename PolygonRange,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(std::ostream& os,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr
#endif
               )
{
  set_ascii_mode(os); // obj is ASCII only
  Generic_writer<std::ostream, File_writer_wavefront> writer(os);
  return writer(points, polygons, np);
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange>
bool write_OBJ(std::ostream& os, const PointRange& points, const PolygonRange& polygons,
               typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return write_OBJ(os, points, polygons, parameters::all_default());
}

/// \endcond

/*!
 * \ingroup PkgStreamSupportIoFuncsOBJ
 *
 * \brief writes the content of `points` and `polygons` in a file named `fname`, using the \ref IOStreamOBJ.
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
 * \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
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
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(const std::string& fname,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr
#endif
               )
{
  std::ofstream os(fname);
  CGAL::IO::set_mode(os, CGAL::IO::ASCII);

  return write_OBJ(os, points, polygons, np);
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange>
bool write_OBJ(const std::string& fname, const PointRange& points, const PolygonRange& polygons,
               typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return write_OBJ(fname, points, polygons, parameters::all_default());
}

/// \endcond

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_OBJ_H
