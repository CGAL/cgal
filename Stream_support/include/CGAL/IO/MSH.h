// Copyright (c) 2025  GeometryFactory Sarl (France).
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
// Author(s)     : Rajdeep Singh

#ifndef CGAL_IO_MSH_H
#define CGAL_IO_MSH_H

#include <CGAL/IO/helpers.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/range/value_type.hpp>
#include <boost/range/size.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace CGAL {

namespace IO {

// ============================================================================
// MSH 2.2 element type codes for surface elements
// ============================================================================

namespace internal {

static const int MSH_LINE         = 1;  // 2-node line (skip)
static const int MSH_TRIANGLE     = 2;  // 3-node triangle
static const int MSH_QUAD         = 3;  // 4-node quadrangle
static const int MSH_POINT        = 15; // 1-node point  (skip)

// Nodes per element type (for the types we care about)
inline int msh_nodes_for_type(int etype)
{
  switch(etype)
  {
    case MSH_LINE:      return 2;
    case MSH_TRIANGLE:  return 3;
    case MSH_QUAD:      return 4;
    case MSH_POINT:     return 1;
    // higher-order triangles / quads
    case 9:             return 6;   // 6-node triangle
    case 10:            return 9;   // 9-node quadrangle
    case 16:            return 8;   // 8-node quadrangle
    case 21:            return 10;  // 10-node triangle
    // volume elements (skip)
    case 4:             return 4;   // tet
    case 5:             return 8;   // hex
    case 6:             return 6;   // prism
    case 7:             return 5;   // pyramid
    case 11:            return 10;  // 10-node tet
    case 17:            return 20;  // 20-node hex
    // default: unknown, we'll skip these
    default:            return -1;
  }
}

} // namespace internal

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

/*!
 * \ingroup PkgStreamSupportIoFuncsMSH
 *
 * \brief reads the content of `is` into `points` and `polygons`, using the
 *        Gmsh \ref IOStreamMSH format (version 2.2, ASCII).
 *
 * Only surface elements (triangles and quadrangles) are extracted; higher-dimensional
 * elements (tetrahedra, hexahedra, …) and lower-dimensional elements (lines, points)
 * are silently skipped.  Mixed meshes that contain both surface and volume elements are
 * therefore handled correctly: only the boundary/surface faces end up in the polygon soup.
 *
 * \attention The polygon soup is not cleared, and the data from the stream are appended.
 *
 * \tparam PointRange a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                    whose value type is the point type
 * \tparam PolygonRange a model of the concepts `SequenceContainer` and `BackInsertionSequence`
 *                      whose `value_type` is itself a model of the concepts `SequenceContainer`
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
template <typename PointRange, typename PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_MSH(std::istream& is,
              PointRange& points,
              PolygonRange& polygons,
              const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
              , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
              )
{
  typedef typename boost::range_value<PointRange>::type    Point;
  typedef typename boost::range_value<PolygonRange>::type  Polygon;

  const bool verbose = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::verbose), false);

  if(!is.good())
  {
    if(verbose)
      std::cerr << "Error: stream is not readable." << std::endl;
    return false;
  }

  // ---- State flags --------------------------------------------------------
  bool found_mesh_format = false;
  bool found_nodes       = false;
  bool found_elements    = false;
  bool binary_mode       = false;

  // gmsh node indices are 1-based.  Map gmsh_tag → 0-based index in `points`.
  std::unordered_map<std::size_t, std::size_t> node_index;

  std::string line;
  while(std::getline(is, line))
  {
    // ---- $MeshFormat -------------------------------------------------------
    if(line == "$MeshFormat")
    {
      std::string version_str;
      int file_type = 0;
      int data_size = 8;

      if(!std::getline(is, line)) break;
      std::istringstream iss(line);
      iss >> version_str >> file_type >> data_size;

      if(file_type != 0)
      {
        if(verbose)
          std::cerr << "Warning: binary MSH files are not supported; only ASCII (file-type 0) is handled." << std::endl;
        binary_mode = true;
        // We still set found_mesh_format so the caller gets a clean 'false'
      }

      // Accept versions 2.x
      double version = 2.0; // default if parsing fails
      try { version = std::stod(version_str); }
      catch(const std::invalid_argument&) {}
      catch(const std::out_of_range&) {}
      if(version < 2.0 || version >= 3.0)
      {
        if(verbose)
          std::cerr << "Warning: MSH version " << version_str
                    << " detected; only version 2.x is supported." << std::endl;
        // still attempt to parse
      }

      found_mesh_format = true;

      // consume $EndMeshFormat
      if(!std::getline(is, line)) break; // "$EndMeshFormat"
      continue;
    }

    // ---- bail early if binary -----------------------------------------------
    if(binary_mode)
      continue;

    // ---- $PhysicalNames (optional, skip) ------------------------------------
    if(line == "$PhysicalNames")
    {
      // skip until $EndPhysicalNames
      while(std::getline(is, line) && line != "$EndPhysicalNames")
        ;
      continue;
    }

    // ---- $Nodes -------------------------------------------------------------
    if(line == "$Nodes")
    {
      std::size_t num_nodes = 0;
      if(!std::getline(is, line)) break;
      {
        std::istringstream iss(line);
        iss >> num_nodes;
      }

      bool nodes_ok = true;
      for(std::size_t i = 0; i < num_nodes; ++i)
      {
        if(!std::getline(is, line))
        {
          if(verbose)
            std::cerr << "Error: unexpected end of file inside $Nodes." << std::endl;
          nodes_ok = false;
          break;
        }
        std::istringstream iss(line);
        std::size_t tag;
        double x, y, z;
        if(!(iss >> tag >> x >> y >> z))
        {
          if(verbose)
            std::cerr << "Error: malformed node line: " << line << std::endl;
          nodes_ok = false;
          break;
        }
        node_index[tag] = points.size();
        Point p;
        internal::fill_point(x, y, z, 1.0, p);
        points.push_back(p);
      }

      if(!nodes_ok)
        return false;

      // consume $EndNodes
      if(!std::getline(is, line)) break; // "$EndNodes"
      found_nodes = true;
      continue;
    }

    // ---- $Elements ----------------------------------------------------------
    if(line == "$Elements")
    {
      std::size_t num_elements = 0;
      if(!std::getline(is, line)) break;
      {
        std::istringstream iss(line);
        iss >> num_elements;
      }

      bool elements_ok = true;
      for(std::size_t i = 0; i < num_elements; ++i)
      {
        if(!std::getline(is, line))
        {
          if(verbose)
            std::cerr << "Error: unexpected end of file inside $Elements." << std::endl;
          elements_ok = false;
          break;
        }

        std::istringstream iss(line);
        std::size_t elem_tag;
        int elem_type, num_tags;

        if(!(iss >> elem_tag >> elem_type >> num_tags))
        {
          if(verbose)
            std::cerr << "Error: malformed element line: " << line << std::endl;
          elements_ok = false;
          break;
        }

        // Skip the tags (physical group, elementary entity, …)
        for(int t = 0; t < num_tags; ++t)
        {
          int tag_val;
          if(!(iss >> tag_val))
          {
            if(verbose)
              std::cerr << "Error: could not read tag " << t
                        << " for element " << elem_tag << "." << std::endl;
            elements_ok = false;
            break;
          }
        }
        if(!elements_ok)
          break;

        int node_count = internal::msh_nodes_for_type(elem_type);
        if(node_count < 0)
        {
          // Unknown element type — we can still try to parse by reading until end of line
          // (already consumed by getline), just skip.
          if(verbose)
            std::cerr << "Warning: unknown element type " << elem_type
                      << " (element " << elem_tag << "), skipping." << std::endl;
          continue;
        }

        // Read node indices for this element
        std::vector<std::size_t> elem_nodes(node_count);
        bool node_ok = true;
        for(int n = 0; n < node_count; ++n)
        {
          std::size_t nid;
          if(!(iss >> nid))
          {
            if(verbose)
              std::cerr << "Error: could not read node indices for element "
                        << elem_tag << "." << std::endl;
            node_ok = false;
            break;
          }

          auto it = node_index.find(nid);
          if(it == node_index.end())
          {
            if(verbose)
              std::cerr << "Error: element " << elem_tag
                        << " references unknown node " << nid << "." << std::endl;
            node_ok = false;
            break;
          }
          elem_nodes[n] = it->second;
        }

        if(!node_ok)
        {
          elements_ok = false;
          break;
        }

        // Only add surface elements to the polygon soup
        if(elem_type == internal::MSH_TRIANGLE || elem_type == internal::MSH_QUAD)
        {
          Polygon face;
          for(int n = 0; n < node_count; ++n)
            face.push_back(static_cast<typename Polygon::value_type>(elem_nodes[n]));
          polygons.push_back(face);
          found_elements = true;
        }
        // else: silently skip lines, points, volume elements, higher-order elements
      }

      if(!elements_ok)
        return false;

      // consume $EndElements
      if(!std::getline(is, line)) break; // "$EndElements"
      continue;
    }

    // ---- Any other $Section — skip until matching $EndXxx ------------------
    if(!line.empty() && line[0] == '$' && line.substr(0, 4) != "$End")
    {
      const std::string end_marker = "$End" + line.substr(1);
      while(std::getline(is, line) && line != end_marker)
        ;
    }
  }

  if(binary_mode)
  {
    if(verbose)
      std::cerr << "Error: binary MSH is not supported." << std::endl;
    return false;
  }

  if(!found_mesh_format)
  {
    if(verbose)
      std::cerr << "Error: $MeshFormat section not found." << std::endl;
    return false;
  }

  if(!found_nodes)
  {
    if(verbose)
      std::cerr << "Error: $Nodes section not found." << std::endl;
    return false;
  }

  if(!found_elements && verbose)
    std::cerr << "Warning: no surface elements (triangles/quads) found." << std::endl;

  return !is.bad();
}

/*!
 * \ingroup PkgStreamSupportIoFuncsMSH
 *
 * \brief reads the content of a file named `fname` into `points` and `polygons`,
 *        using the Gmsh \ref IOStreamMSH format (version 2.2, ASCII).
 *
 * \attention The polygon soup is not cleared, and the data from the file are appended.
 *
 * \tparam PointRange a model of the concepts `RandomAccessContainer` and `BackInsertionSequence`
 *                    whose value type is the point type
 * \tparam PolygonRange a model of the concepts `SequenceContainer` and `BackInsertionSequence`
 *                      whose `value_type` is itself a model of the concepts `SequenceContainer`
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
template <typename PointRange, typename PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_MSH(const std::string& fname,
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

  if(!is)
  {
    const bool verbose = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::verbose), false);
    if(verbose)
      std::cerr << "Error: cannot open file '" << fname << "'." << std::endl;
    return false;
  }

  return read_MSH(is, points, polygons, np);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 * \ingroup PkgStreamSupportIoFuncsMSH
 *
 * \brief writes the content of `points` and `polygons` in `os`,
 *        using the Gmsh \ref IOStreamMSH format (version 2.2, ASCII).
 *
 * Triangles are written as MSH element type 2 and quadrangles as element type 3.
 * Polygons with more than four vertices are fan-triangulated before writing.
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
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
 *     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{the precision of the stream `os`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if the writing was successful, `false` otherwise.
 */
template <typename PointRange, typename PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_MSH(std::ostream& os,
               const PointRange& points,
               const PolygonRange& polygons,
               const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<internal::is_Range<PolygonRange>::value>* = nullptr
#endif
               )
{
  typedef typename boost::range_value<PointRange>::type   Point;
  typedef typename boost::range_value<PolygonRange>::type Polygon;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, CGAL_NP_CLASS>::type PointMap;
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));

  if(!os.good())
    return false;

  set_stream_precision_from_NP(os, np);

  // ---- $MeshFormat --------------------------------------------------------
  os << "$MeshFormat\n";
  os << "2.2 0 8\n";
  os << "$EndMeshFormat\n";

  // ---- $Nodes -------------------------------------------------------------
  const std::size_t num_nodes = boost::size(points);
  os << "$Nodes\n";
  os << num_nodes << "\n";
  {
    std::size_t tag = 1; // gmsh is 1-indexed
    for(std::size_t i = 0; i < num_nodes; ++i, ++tag)
    {
      const Point& p = get(point_map, *(std::begin(points) + i));
      os << tag << " "
         << CGAL::to_double(p.x()) << " "
         << CGAL::to_double(p.y()) << " "
         << CGAL::to_double(p.z()) << "\n";
    }
  }
  os << "$EndNodes\n";

  // ---- $Elements ----------------------------------------------------------
  // Pre-compute the number of output elements.
  // Quads stay as quads; polygons with > 4 vertices are fan-triangulated.
  std::size_t num_output_elements = 0;
  for(const Polygon& face : polygons)
  {
    const std::size_t n = face.size();
    if(n == 3)       num_output_elements += 1;           // one triangle
    else if(n == 4)  num_output_elements += 1;           // one quad
    else if(n > 4)   num_output_elements += (n - 2);     // fan triangulation
    // n < 3: degenerate, skip
  }

  os << "$Elements\n";
  os << num_output_elements << "\n";

  std::size_t elem_tag = 1;
  for(const Polygon& face : polygons)
  {
    const std::size_t n = face.size();
    if(n < 3)
      continue; // skip degenerate faces silently

    if(n == 3)
    {
      // MSH type 2: triangle, 2 default tags (physical=0, elementary=0)
      os << elem_tag++ << " 2 2 0 0 "
         << (face[0] + 1) << " " << (face[1] + 1) << " " << (face[2] + 1) << "\n";
    }
    else if(n == 4)
    {
      // MSH type 3: quadrangle
      os << elem_tag++ << " 3 2 0 0 "
         << (face[0] + 1) << " " << (face[1] + 1) << " "
         << (face[2] + 1) << " " << (face[3] + 1) << "\n";
    }
    else
    {
      // Fan triangulation: all triangles share face[0]
      for(std::size_t k = 1; k + 1 < n; ++k)
      {
        os << elem_tag++ << " 2 2 0 0 "
           << (face[0] + 1) << " " << (face[k] + 1) << " " << (face[k + 1] + 1) << "\n";
      }
    }
  }
  os << "$EndElements\n";

  return !os.fail();
}

/*!
 * \ingroup PkgStreamSupportIoFuncsMSH
 *
 * \brief writes the content of `points` and `polygons` in a file named `fname`,
 *        using the Gmsh \ref IOStreamMSH format (version 2.2, ASCII).
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
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
 *     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
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
template <typename PointRange, typename PolygonRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_MSH(const std::string& fname,
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

  if(!os)
  {
    const bool verbose = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::verbose), false);
    if(verbose)
      std::cerr << "Error: cannot open file '" << fname << "' for writing." << std::endl;
    return false;
  }

  return write_MSH(os, points, polygons, np);
}

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_MSH_H
