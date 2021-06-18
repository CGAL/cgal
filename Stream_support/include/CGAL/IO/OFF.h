// Copyright (c) 2015-2020 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau and Sebastien Loriot

#ifndef CGAL_IO_OFF_H
#define CGAL_IO_OFF_H

#include <CGAL/IO/OFF/Scanner_OFF.h>
#include <CGAL/IO/OFF/File_scanner_OFF.h>
#include <CGAL/IO/OFF/File_writer_OFF.h>
#include <CGAL/IO/OFF/generic_copy_OFF.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/IO/Generic_writer.h>

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/use.h>

#include <boost/range/value_type.hpp>
#include <boost/utility/enable_if.hpp>

#include <fstream>
#include <iostream>
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

//help fixing conversion warnings without using value_type, which is forbidden by c++20
template< typename S, typename T>
void integer_type_converter(S& p1, const T& p2)
{
  p1 = static_cast<S>(p2);
}

template <typename PointRange, typename PolygonRange,
          typename VertexNormalOutputIterator,
          typename VertexColorOutputIterator,
          typename VertexTextureOutputIterator,
          typename FaceColorOutputIterator>
bool read_OFF(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              VertexNormalOutputIterator vn_out,
              VertexColorOutputIterator vc_out,
              VertexTextureOutputIterator vt_out,
              FaceColorOutputIterator fc_out,
              const bool verbose = false)
{
  typedef typename boost::range_value<PointRange>::type                               Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel                                 Kernel;
  typedef typename Kernel::Point_2                                                    Texture;
  typedef typename Kernel::Vector_3                                                   Normal;
  typedef typename Kernel::FT                                                         FT;
  typedef CGAL::IO::Color                                                                 Color;


  if(!is.good()){
    if(verbose)
      std::cerr<<"File doesn't exist."<<std::endl;
    return false;
  }

  CGAL::File_scanner_OFF scanner(is);
  if(is.fail())
    return false;
  points.resize(scanner.size_of_vertices());
  polygons.resize(scanner.size_of_facets());

  for(std::size_t i=0; i<scanner.size_of_vertices(); ++i)
  {
    double x(0), y(0), z(0), w(0);
    scanner.scan_vertex(x, y, z, w);
    CGAL_assertion(w != 0);
    internal::fill_point(x, y, z, w, points[i]);

    if(scanner.has_normals())
    {
      double nx, ny, nz;
      scanner.scan_normal(nx, ny, nz);
      *vn_out++ = Normal(FT(nx), FT(ny), FT(nz));
    }

    if(scanner.has_vcolors())
    {
      unsigned char r=0, g=0, b=0;
      scanner.scan_color(r, g, b);
      *vc_out++ = Color(r,g,b);
    }

    if(scanner.has_textures())
    {
      double nx, ny;
      scanner.scan_texture(nx, ny);
      *vt_out++ = Texture(FT(nx), FT(ny));
    }
    if(!is.good())
      return false;
  }


  for(std::size_t i=0; i<scanner.size_of_facets(); ++i)
  {
    std::size_t no(-1);
    scanner.scan_facet(no, i);

    if((!is.eof() && !is.good()) || no == std::size_t(-1))
      return false;

    CGAL::internal::resize(polygons[i], no);
    for(std::size_t j=0; j<no; ++j)
    {
      std::size_t id = 0;
      scanner.scan_facet_vertex_index(id,j+1, i);
      if(!is)
      {
        return false;
      }
      if(id < scanner.size_of_vertices())
        integer_type_converter(polygons[i][j], id);
      else
        return false;
    }
    if(scanner.has_fcolors())
    {
      unsigned char r=0, g=0, b=0;
      scanner.scan_color(r,g,b);
      *fc_out++ = Color(r,g,b);
    }
  }

  return !is.fail();
}

} // namespace internal

/*!
 * \ingroup PkgStreamSupportIoFuncsOFF
 *
 * \brief reads the content of `is` into `points` and `polygons`, using the \ref IOStreamOFF.
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
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{indicates whether output warnings and error messages should be printed or not.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the reading was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
              , typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr
#endif
              )
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  return internal::read_OFF(is, points, polygons,
                            choose_parameter(get_parameter(np, internal_np::vertex_normal_output_iterator),
                                             CGAL::Emptyset_iterator()),
                            choose_parameter(get_parameter(np, internal_np::vertex_color_output_iterator),
                                             CGAL::Emptyset_iterator()),
                            choose_parameter(get_parameter(np, internal_np::vertex_texture_output_iterator),
                                             CGAL::Emptyset_iterator()),
                            choose_parameter(get_parameter(np, internal_np::face_color_output_iterator),
                                             CGAL::Emptyset_iterator()),
                            choose_parameter(get_parameter(np, internal_np::verbose), true));
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange>
bool read_OFF(std::istream& is, PointRange& points, PolygonRange& polygons,
              typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return read_OFF(is, points, polygons, parameters::all_default());
}

/// \endcond

/*!
 * \ingroup PkgStreamSupportIoFuncsOFF
 *
 * \brief reads the content of the file `fname` into `points` and `polygons`, using the \ref IOStreamOFF.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 * \tparam PolygonRange a model of the concepts `SequenceContainer` and `BackInsertionSequence`
 *                      whose `value_type` is itself a model of the concept `SequenceContainer`
 *                      and `BackInsertionSequence` whose `value_type` is an unsigned integer type
 *                      convertible to `std::size_t`
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
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the reading was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const std::string& fname,
              PointRange& points,
              PolygonRange& polygons,
              const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
              , typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr
#endif
              )
{
  std::ifstream in(fname);
  return read_OFF(in, points, polygons, np);
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange>
bool read_OFF(const std::string& fname, PointRange& points, PolygonRange& polygons,
              typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return read_OFF(fname, points, polygons, parameters::all_default());
}

/// \endcond

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 * \ingroup PkgStreamSupportIoFuncsOFF
 *
 * \brief writes the content of `points` and `polygons` in `os`, using the \ref IOStreamOFF.
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
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(std::ostream& os,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr
#endif
               )
{
  Generic_writer<std::ostream, File_writer_OFF> writer(os);
  return writer(points, polygons, np);
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange>
bool write_OFF(std::ostream& os, const PointRange& points, const PolygonRange& polygons
               , typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return write_OFF(os, points, polygons, parameters::all_default());
}

/// \endcond

/*!
 * \ingroup PkgStreamSupportIoFuncsOFF
 *
 * \brief writes the content of `points` and `polygons` in the file `fname`, using the \ref IOStreamOFF.
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
template <typename PointRange, typename PolygonRange, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(const std::string& fname,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr
#endif
               )
{
  std::ofstream os(fname);
  Generic_writer<std::ostream, File_writer_OFF> writer(os);
  return writer(points, polygons, np);
}

/// \cond SKIP_IN_MANUAL

template <typename PointRange, typename PolygonRange>
bool write_OFF(const std::string& fname, const PointRange& points, const PolygonRange& polygons,
               typename boost::enable_if<internal::is_Range<PolygonRange> >::type* = nullptr)
{
  return write_OFF(fname, points, polygons, parameters::all_default());
}

/// \endcond

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_OFF_H
