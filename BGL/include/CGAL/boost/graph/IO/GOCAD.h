// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//                 Mael Rouxel-Labbé

#ifndef CGAL_BGL_IO_GOCAD_H
#define CGAL_BGL_IO_GOCAD_H

#include <CGAL/IO/GOCAD.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/IO/Generic_facegraph_builder.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

#include <boost/container/flat_map.hpp>

#include <fstream>
#include <iostream>
#include <utility>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#endif

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

namespace IO {
namespace internal {

// Use CRTP to gain access to the protected members without getters/setters.
template <typename FaceGraph, typename Point>
class GOCAD_builder
  : public Generic_facegraph_builder<FaceGraph, Point, GOCAD_builder<FaceGraph, Point> >
{
  typedef GOCAD_builder<FaceGraph, Point>                                       Self;
  typedef Generic_facegraph_builder<FaceGraph, Point, Self>                     Base;

  typedef typename Base::Point_container                                        Point_container;
  typedef typename Base::Face                                                   Face;
  typedef typename Base::Face_container                                         Face_container;

public:
  GOCAD_builder(std::istream& is, bool verbose) : Base(is, verbose) { }

  // @check ascii
  template <typename NamedParameters>
  bool read(std::istream& is,
            Point_container& points,
            Face_container& faces,
            const NamedParameters& np,
            bool verbose)
  {
    std::pair<std::string, std::string> name_and_color;
    bool res = read_GOCAD(is, name_and_color, points, faces, np, verbose);
    if(res)
    {
      name = name_and_color.first;
      color = name_and_color.second;
    }
    return res;
  }

public:
  std::string name;
  std::string color;
};

} // namespace internal
} // namespace IO

/// \ingroup PkgBGLIOFct
///
/// reads the graph `g` from the input stream in the TS format.
///
/// \tparam FaceGraph a model of `MutableFaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param is the input stream
/// \param name_and_color name and color of the mesh
/// \param g the graph to be built from the input data
/// \param verbose whether extra information is printed when an incident occurs during reading
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `g`}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `FaceGraph`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \pre The data must represent a 2-manifold
///
/// \attention The graph `g` is not cleared, and the data from the stream is added.
///
/// \returns `true` if the resulting mesh is valid.
///
/// \see \ref IOStreamGocad
template <typename FaceGraph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(std::istream& is,
                std::pair<std::string, std::string>& name_and_color,
                FaceGraph& g,
                const CGAL_BGL_NP_CLASS& np,
                bool verbose = true)
{
  typedef typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::type VPM;
  typedef typename boost::property_traits<VPM>::value_type                     Point;

  IO::internal::GOCAD_builder<FaceGraph, Point> builder(is, verbose);
  if(!builder(g, np))
    return false;

  name_and_color.first = builder.name;
  name_and_color.second = builder.color;

  return is_valid(g); // @fixme keep validity check?
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(std::istream& is, FaceGraph& g, const CGAL_BGL_NP_CLASS& np, bool verbose = true)
{
  std::pair<std::string, std::string> dummy;
  return read_GOCAD(is, dummy, g, np, verbose);
}

template <typename FaceGraph>
bool read_GOCAD(std::istream& is, std::pair<std::string, std::string>& name_and_color, FaceGraph& g)
{
  return read_GOCAD(is, name_and_color, g, parameters::all_default());
}

/// \ingroup PkgBGLIOFct
///
/// reads the graph `g` from the file `fname` in the TS format.
///
/// \tparam FaceGraph a model of `MutableFaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param fname the name of the input file
/// \param name_and_color name and color of the mesh
/// \param g the graph to be built from the input data
/// \param verbose whether extra information is printed when an incident occurs during reading
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `g`}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `FaceGraph`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \sa Overloads of this function for specific models of the concept `FaceGraph`.
///
/// \pre The data must represent a 2-manifold
///
/// \see \ref IOStreamGocad
template <typename FaceGraph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(const char* fname,
                std::pair<std::string, std::string>& name_and_color,
                FaceGraph& g,
                const CGAL_BGL_NP_CLASS& np,
                bool verbose = true)
{
  std::ifstream is(fname);
  return read_GOCAD(is, name_and_color, g, np, verbose);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(const char* fname, FaceGraph& g, const CGAL_BGL_NP_CLASS& np, bool verbose = true)
{
  std::pair<std::string, std::string> dummy;
  return read_GOCAD(fname, dummy, g, np, verbose);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(const std::string& fname, FaceGraph& g, CGAL_BGL_NP_CLASS np, bool verbose = true)
{
  return read_GOCAD(fname.c_str(), g, np, verbose);
}

template <typename FaceGraph>
bool read_GOCAD(std::istream& is, FaceGraph& g) { return read_GOCAD(is, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_GOCAD(const char* fname, FaceGraph& g) { return read_GOCAD(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_GOCAD(const std::string& fname, FaceGraph& g) { return read_GOCAD(fname, g, parameters::all_default()); }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/// \ingroup PkgBGLIOFct
///
/// writes the graph `g` in the TS format into `os`.
///
/// \tparam FaceGraph a model of `FaceListGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param os the stream into which `g` is dumped
/// \param name the name that will be assigned to `g` in the file
/// \param g the graph to be output
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `g`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `FaceGraph`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \see \ref IOStreamGocad
template <typename FaceGraph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(std::ostream& os,
                 const char* name,
                 const FaceGraph& g,
                 const CGAL_BGL_NP_CLASS& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor      face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type   vertices_size_type;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  if(!os.good())
    return false;

  os << "GOCAD TSurf 1\n"
        "HEADER {\n"
        "name:";
  os << name << std::endl;
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

  boost::container::flat_map<vertex_descriptor, vertices_size_type> reindex;

  vertices_size_type i = 0;
  for(const vertex_descriptor v : vertices(g))
  {
    os << "VRTX " << i << " " << get(vpm, v) << "\n";
    reindex[v] = i++;
  }

  for(const face_descriptor f : faces(g))
  {
    halfedge_descriptor h = halfedge(f, g);
    os << "TRGL " << reindex[target(prev(h, g), g)] << " "
                  << reindex[target(h, g)] << " "
                  << reindex[target(next(h, g), g)] << "\n";
  }

  os << "END" << std::endl;

  return os.good();
}

template <typename FaceGraph>
bool write_GOCAD(std::ostream& os, const char* name, const FaceGraph& g)
{
  return write_GOCAD(os, name, g, parameters::all_default());
}

/// \ingroup PkgBGLIOFct
///
/// writes the graph `g` in the TS format into a file named `fname`. In this overload,
/// `fname` is used as the name of the graph within the file.
///
/// \tparam FaceGraph a model of `FaceListGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param fname the name of the output file
/// \param g the graph to be output
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `g`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `FaceGraph`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \sa Overloads of this function for specific models of the concept `FaceGraph`.
///
/// \see \ref IOStreamGocad
template <typename FaceGraph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(const char* fname,
                 const FaceGraph& g,
                 const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream os(fname);
  return write_GOCAD(os, fname, g, np);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(const std::string& fname, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  return write_GOCAD(fname.c_str(), g, np);
}

/// \ingroup PkgBGLIOFct
///
/// writes the graph `g` in the TS format into `os`. The name
/// assigned to `g`in the output is `anonymous`.
///
/// \tparam FaceGraph a model of `FaceListGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param os the output stream
/// \param g the graph to be output
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `g`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `FaceGraph`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \see \ref IOStreamGocad
template <typename FaceGraph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(std::ostream& os,
                 const FaceGraph& g,
                 const CGAL_BGL_NP_CLASS& np)
{
  return write_GOCAD(os, "anonymous", g, np);
}
template <typename FaceGraph>
bool write_GOCAD(std::ostream& os, const FaceGraph& g)
{
  return write_GOCAD(os, "anonymous", g, parameters::all_default());
}

template <typename FaceGraph>
bool write_GOCAD(const char* fname, const FaceGraph& g) { return write_GOCAD(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_GOCAD(const std::string& fname, const FaceGraph& g) { return write_GOCAD(fname, g, parameters::all_default()); }

} // namespace CGAL

#endif // CGAL_BGL_IO_GOCAD_H
