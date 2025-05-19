// Copyright (c) 2018-2020  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_STL_H
#define CGAL_IO_STL_H

#include <CGAL/IO/STL/STL_reader.h>
#include <CGAL/IO/helpers.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <boost/range/value_type.hpp>

#include <type_traits>
#include <iostream>
#include <fstream>
#include <string>

namespace CGAL {

namespace IO {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

/*!
 * \ingroup PkgStreamSupportIoFuncsSTL
 *
 * \brief reads the content of `is` into `points` and `facets`, using the \ref IOStreamSTL.
 *
 * \attention The polygon soup is not cleared, and the data from the stream are appended.
 *
 * \attention To read a binary file, the flag `std::ios::binary` must be set during the creation of the `ifstream`.
 *
 * \tparam PointRange a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                    whose value type is the point type
 * \tparam TriangleRange a model of the concept `SequenceContainer`
 *                      whose `value_type` is itself a model of the concept `SequenceContainer`
 *                      whose `value_type` is an unsigned integer type convertible to `std::size_t`.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param is the input stream
 * \param points points of the soup of triangles
 * \param facets a range of triangles; each triangle uses the indices of the points in `points`.
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
template <typename PointRange, typename TriangleRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_STL(std::istream& is,
              PointRange& points,
              TriangleRange& facets,
              const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
              , std::enable_if_t<internal::is_Range<TriangleRange>::value>* = nullptr
#endif
              )
{
  const bool verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);

  if(!is.good())
  {
    if(verbose)
      std::cerr << "File doesn't exist." << std::endl;
    return false;
  }

  // Ignore all initial whitespace
  unsigned char c;

  while(is.read(reinterpret_cast<char*>(&c), sizeof(c)))
  {
    if(!isspace(c))
    {
      is.unget(); // move back to the first interesting char
      break;
    }
  }

  if(!is.good()) // reached the end
    return true;

  // Read the 5 first characters to check if the first word is "solid"
  std::string s;

  char word[6];
  if(is.read(reinterpret_cast<char*>(&word[0]), sizeof(c)) &&
     is.read(reinterpret_cast<char*>(&word[1]), sizeof(c)) &&
     is.read(reinterpret_cast<char*>(&word[2]), sizeof(c)) &&
     is.read(reinterpret_cast<char*>(&word[3]), sizeof(c)) &&
     is.read(reinterpret_cast<char*>(&word[4]), sizeof(c)) &&
     is.read(reinterpret_cast<char*>(&word[5]), sizeof(c)))
  {
    s = std::string(word, 5);
  }
  else
  {
    return true; // empty file
  }

  // If the first word is not 'solid', the file must be binary
  if(s != "solid" || (word[5] !='\n' && word[5] !='\r' && word[5] != ' '))
  {
    if(internal::parse_binary_STL(is, points, facets, verbose))
    {
      return true;
    }
    else
    {
      // If we failed to read it as a binary, try as ASCII just in case...
      // The file does not start with 'solid' anyway, so it's fine to reset it.
      is.clear();
      is.seekg(0, std::ios::beg);
      return internal::parse_ASCII_STL(is, points, facets, verbose);
    }

  }
  // Now, we have found the keyword "solid" which is supposed to indicate that the file is ASCII
  is.clear();
  is.seekg(0, std::ios::beg); // the parser needs to read all "solid" to work correctly.
  if(internal::parse_ASCII_STL(is, points, facets, verbose))
  {
    // correctly read the input as an ASCII file
    return true;
  }
  else// Failed to read the ASCII file
  {
    // It might have actually have been a binary file... ?
    return internal::parse_binary_STL(is, points, facets, verbose);
  }
}

/*!
 * \ingroup PkgStreamSupportIoFuncsSTL
 *
 * \brief reads the content of a file named `fname` into `points` and `facets`, using the \ref IOStreamSTL.
 *
 *  If `use_binary_mode` is `true`, but the reading fails, \ascii reading will be automatically tested.
 * \attention The polygon soup is not cleared, and the data from the file are appended.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * \tparam TriangleRange a model of the concept `SequenceContainer`
 *                      whose `value_type` is itself a model of the concept `SequenceContainer`
 *                      whose `value_type` is an unsigned integer type convertible to `std::size_t`.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the path to the input file
 * \param points points of the soup of triangles
 * \param facets a range of triangles; each triangle uses the indices of the points in `points`.
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
template <typename PointRange, typename TriangleRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_STL(const std::string& fname,
              PointRange& points,
              TriangleRange& facets,
              const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
              , std::enable_if_t<internal::is_Range<TriangleRange>::value>* = nullptr
#endif
              )
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  const bool binary = parameters::choose_parameter(parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ifstream is(fname, std::ios::binary);
    CGAL::IO::set_mode(is, BINARY);
    if(read_STL(is, points, facets, np))
    {
      return true;
    }
    points.clear();
    facets.clear();
  }
  std::ifstream is(fname);
  CGAL::IO::set_mode(is, CGAL::IO::ASCII);
  bool v = choose_parameter(get_parameter(np, internal_np::verbose),
                            false);
  return read_STL(is, points, facets, CGAL::parameters::verbose(v).use_binary_mode(false));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace internal {

template <typename K>
typename K::Vector_3 construct_normal_of_STL_face(const typename K::Point_3& p,
                                                  const typename K::Point_3& q,
                                                  const typename K::Point_3& r,
                                                  const K& k)
{
  typedef typename K::FT FT;
  typedef typename K::Vector_3 Vector;

  if(k.collinear_3_object()(p, q, r))
    return k.construct_vector_3_object()(1, 0, 0);

  Vector res = k.construct_orthogonal_vector_3_object()(p, q, r);
  const FT sql = k.compute_squared_length_3_object()(res);
  res = k.construct_divided_vector_3_object()(res, CGAL::approximate_sqrt(sql));

  return res;
}

} // namespace internal

/*!
 * \ingroup PkgStreamSupportIoFuncsSTL
 *
 * \brief writes the content of `points` and `facets` in `os`, using the \ref IOStreamSTL.
 *
 * \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation
 *            of the `ofstream`, and the \link PkgStreamSupportEnumRef `IO::Mode` \endlink
 *            of the stream must be set to `BINARY`.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * \tparam TriangleRange a model of the concept `SequenceContainer`
 *                      whose `value_type` is itself a model of the concept `SequenceContainer`
 *                      whose `value_type` is an unsigned integer type convertible to `std::size_t`.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param os the output stream
 * \param points points of the soup of triangles
 * \param facets a range of triangles; each triangle uses the indices of the points in `points`.
 * \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
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
template <typename PointRange, typename TriangleRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_STL(std::ostream& os,
               const PointRange& points,
               const TriangleRange& facets,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<TriangleRange>::value>* = nullptr
#endif
               )
{
  typedef typename boost::range_value<TriangleRange>::type                  Triangle;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, CGAL_NP_CLASS>::type   PointMap;
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));

  typedef typename boost::property_traits<PointMap>::value_type             Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel                       K;
  typedef typename K::Vector_3                                              Vector_3;

  K k = choose_parameter<K>(get_parameter(np, internal_np::geom_traits));

  if(!os.good())
    return false;

  set_stream_precision_from_NP(os, np);

  if(get_mode(os) == BINARY)
  {
    os << "FileType: Binary                                                                ";
    const std::uint32_t N32 = static_cast<std::uint32_t>(facets.size());
    os.write(reinterpret_cast<const char *>(&N32), sizeof(N32));

    for(const Triangle& face : facets)
    {
      const Point& p = get(point_map, points[face[0]]);
      const Point& q = get(point_map, points[face[1]]);
      const Point& r = get(point_map, points[face[2]]);

      const Vector_3 n = collinear(p,q,r) ? Vector_3(1,0,0) : unit_normal(p,q,r);

      const float coords[12] = { static_cast<float>(n.x()), static_cast<float>(n.y()), static_cast<float>(n.z()),
                                 static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()),
                                 static_cast<float>(q.x()), static_cast<float>(q.y()), static_cast<float>(q.z()),
                                 static_cast<float>(r.x()), static_cast<float>(r.y()), static_cast<float>(r.z()) };

      for(int i=0; i<12; ++i)
        os.write(reinterpret_cast<const char *>(&coords[i]), sizeof(coords[i]));
      os << "  ";
    }
  }
  else
  {
    os << "solid\n";
    for(const Triangle& face : facets)
    {
      const Point& p = get(point_map, points[face[0]]);
      const Point& q = get(point_map, points[face[1]]);
      const Point& r = get(point_map, points[face[2]]);

      const Vector_3 n = internal::construct_normal_of_STL_face(p, q, r, k);
      os << "facet normal " << n << "\nouter loop\n";
      os << "vertex " << p << "\n";
      os << "vertex " << q << "\n";
      os << "vertex " << r << "\n";
      os << "endloop\nendfacet\n";
    }
    os << "endsolid"<<std::endl;
  }

  return !os.fail();
}

/*!
 * \ingroup PkgStreamSupportIoFuncsSTL
 *
 * \brief writes the content of `points` and `facets` in a file named `fname`, using the \ref IOStreamSTL.
 *
 * \attention The polygon soup is not cleared, and the data from the file are appended.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * \tparam TriangleRange a model of the concept `SequenceContainer`
 *                      whose `value_type` is itself a model of the concept `SequenceContainer`
 *                      whose `value_type` is an unsigned integer type convertible to `std::size_t`.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the path to the output file
 * \param points points of the soup of triangles
 * \param facets a range of triangles; each triangle uses the indices of the points in `points`.
 * \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
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
template <typename PointRange, typename TriangleRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_STL(const std::string& fname,
               const PointRange& points,
               const TriangleRange& facets,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<TriangleRange>::value>* = nullptr
#endif
               )
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ofstream os(fname, std::ios::binary);
    CGAL::IO::set_mode(os, CGAL::IO::BINARY);
    return write_STL(os, points, facets, np);
  }
  else
  {
    std::ofstream os(fname);
    CGAL::IO::set_mode(os, CGAL::IO::ASCII);
    return write_STL(os, points, facets, np);
  }
}

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_STL_H
