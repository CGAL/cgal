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

#ifndef CGAL_BGL_IO_PLY_H
#define CGAL_BGL_IO_PLY_H

#include <CGAL/IO/PLY.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/boost/graph/IO/Generic_facegraph_builder.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/utility/enable_if.hpp>

#include <fstream>
#include <string>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_typename NamedParameters
#endif

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

namespace IO {
namespace internal {

// Use CRTP to gain access to the protected members without getters/setters.
template <typename Graph, typename Point>
class PLY_builder
  : public Generic_facegraph_builder<Graph, Point, PLY_builder<Graph, Point> >
{
  typedef PLY_builder<Graph, Point>                                         Self;
  typedef Generic_facegraph_builder<Graph, Point, Self>                     Base;

  typedef typename Base::Point_container                                    Point_container;
  typedef typename Base::Face                                               Face;
  typedef typename Base::Face_container                                     Face_container;

public:
  PLY_builder(std::istream& is, bool verbose) : Base(is, verbose) { }

  template <typename NamedParameters>
  bool read(std::istream& is,
            Point_container& points,
            Face_container& faces,
            const NamedParameters& np,
            bool verbose)
  {
    return read_PLY(is, points, faces, np, verbose);
  }
};

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_PLY_BGL(std::istream& is,
                  Graph& g,
                  const CGAL_BGL_NP_CLASS& np,
                  bool verbose = true)
{
  typedef typename CGAL::GetVertexPointMap<Graph, CGAL_BGL_NP_CLASS>::type      VPM;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  IO::internal::PLY_builder<Graph, Point> builder(is, verbose);
  return builder(g, np);
}

} // namespace internal
} // namespace IO

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the PLY format.

  \tparam Graph a model of `MutableFaceGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param is the input stream
  \param g the graph to be built from the input data
  \param verbose whether extra information is printed when an incident occurs during reading
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd

   \cgalParamNBegin{vertex_index_map}
     \cgalParamDescription{a property map associating to each vertex of `graph` a unique index}
     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                    as key type and `std::size_t` as value type}
     \cgalParamDefault{vertex indices that may exist in the input will be ignored}
   \cgalParamNEnd

    \cgalParamNBegin{vertex_color_map}
      \cgalParamDescription{a property map associating colors to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{vertex colors that may exist in the input will be ignored}
    \cgalParamNEnd

    \cgalParamNBegin{face_color_map}
      \cgalParamDescription{a property map associating colors to the faces of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%face_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{face colors that may exist in the input will be ignored}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamPLY
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(std::istream& is,
              Graph& g,
              const CGAL_BGL_NP_CLASS& np,
              bool verbose = true
#ifndef DOXYGEN_RUNNING
              , typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
              )
{
  return IO::internal::read_PLY_BGL(is, g, np, verbose);
}

template <typename Graph>
bool read_PLY(std::istream& is, Graph& g,
              typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return IO::internal::read_PLY_BGL(is, g, parameters::all_default());
}

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from a file named `fname`, in the PLY format.

  \tparam Graph a model of `MutableFaceGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the name of the input file
  \param g the graph to be built from the input data
  \param verbose whether extra information is printed when an incident occurs during reading
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd

   \cgalParamNBegin{vertex_index_map}
     \cgalParamDescription{a property map associating to each vertex of `graph` a unique index}
     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                    as key type and `std::size_t` as value type}
     \cgalParamDefault{vertex indices that may exist in the input will be ignored}
   \cgalParamNEnd

    \cgalParamNBegin{vertex_color_map}
      \cgalParamDescription{a property map associating colors to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{vertex colors that may exist in the input will be ignored}
    \cgalParamNEnd

    \cgalParamNBegin{face_color_map}
      \cgalParamDescription{a property map associating colors to the faces of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%face_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{face colors that may exist in the input will be ignored}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamPLY
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const char* fname,
              Graph& g,
              const CGAL_BGL_NP_CLASS& np,
              bool verbose = true
#ifndef DOXYGEN_RUNNING
              , typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
              )
{
  std::ifstream is(fname);
  return IO::internal::read_PLY_BGL(is, g, np, verbose);
}

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const std::string& fname, Graph& g, const CGAL_BGL_NP_CLASS& np, bool verbose = true,
              typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  std::ifstream is(fname.c_str());
  return IO::internal::read_PLY_BGL(is, g, np, verbose);
}

template <typename Graph>
bool read_PLY(const char* fname, Graph& g,
              typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  std::ifstream is(fname);
  return IO::internal::read_PLY_BGL(is, g, parameters::all_default());
}

template <typename Graph>
bool read_PLY(const std::string& fname, Graph& g,
              typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return IO::internal::read_PLY_BGL(fname.c_str(), g, parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
 \ingroup PkgBGLIOFct

 writes the graph in an output stream in PLY format.

 \tparam Graph a model of `FaceListGraph`
 \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

 \param os the output stream
 \param g the graph to be output
 \param comments a string included line by line in the header of the PLY stream (each line will be precedeed by "comment ")
 \param np optional \ref bgl_namedparameters "Named Parameters" described below

 \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd

   \cgalParamNBegin{vertex_index_map}
     \cgalParamDescription{a property map associating to each vertex of `graph` a unique index}
     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                    as key type and `std::size_t` as value type}
     \cgalParamDefault{vertex indices that may exist in the input will be ignored}
   \cgalParamNEnd

    \cgalParamNBegin{vertex_color_map}
      \cgalParamDescription{a property map associating colors to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{vertex colors that may exist in the input will be ignored}
    \cgalParamNEnd

    \cgalParamNBegin{face_color_map}
      \cgalParamDescription{a property map associating colors to the faces of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%face_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{face colors that may exist in the input will be ignored}
    \cgalParamNEnd
 \cgalNamedParamsEnd

 \returns `true` if writing was successful.

 \see \ref IOStreamPLY
*/
template <typename Graph, typename NamedParameters>
bool write_PLY(std::ostream& os,
               const Graph& g,
               const std::string& comments,
               const NamedParameters& np
#ifndef DOXYGEN_RUNNING
               , typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
               )
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor                          vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor                        halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor                            face_descriptor;

  typedef typename CGAL::GetInitializedVertexIndexMap<Graph, NamedParameters>::const_type VIMap;
  typedef typename GetVertexPointMap<Graph, NamedParameters>::const_type                  Vpm;
  typedef typename boost::property_traits<Vpm>::reference                                 Point_3;
  typedef CGAL::Color                                                                     Color;
  typedef typename internal_np::Lookup_named_param_def<
                     internal_np::vertex_color_map_t,
                     NamedParameters,
                     Constant_property_map<vertex_descriptor, Color> >::type              VCM;
  typedef typename internal_np::Lookup_named_param_def<
                     internal_np::face_color_map_t,
                     NamedParameters,
                     Constant_property_map<face_descriptor, Color> >::type                FCM;

  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;

  VCM vcm = choose_parameter(get_parameter(np, internal_np::vertex_color_map), VCM());
  FCM fcm = choose_parameter(get_parameter(np, internal_np::face_color_map), FCM());

  bool has_vcolor = !is_default_parameter(get_parameter(np, internal_np::vertex_color_map));
  bool has_fcolor = !is_default_parameter(get_parameter(np, internal_np::face_color_map));
  VIMap vim = CGAL::get_initialized_vertex_index_map(g, np);
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(boost::vertex_point, g));

  if(!os.good())
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // Write header
  os << "ply" << std::endl
      << ((get_mode(os) == IO::BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
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

  os << "element vertex " << num_vertices(g) << std::endl;
  IO::internal::output_property_header(os, make_ply_point_writer (CGAL::Identity_property_map<Point_3>()));
  //if vcm is not default add v:color property
  if(has_vcolor)
  {
    os << "property uchar red" << std::endl
       << "property uchar green" << std::endl
       << "property uchar blue" << std::endl
       << "property uchar alpha" << std::endl;
  }

  os << "element face " << num_faces(g) << std::endl;
  IO::internal::output_property_header(
        os, std::make_pair(CGAL::Identity_property_map<std::vector<std::size_t> >(),
                            PLY_property<std::vector<int> >("vertex_indices")));
  //if fcm is not default add f:color property
  if(has_fcolor)
  {
    os << "property uchar red" << std::endl
       << "property uchar green" << std::endl
       << "property uchar blue" << std::endl
       << "property uchar alpha" << std::endl;
  }
  os << "end_header" << std::endl;

  for(vertex_descriptor vd : vertices(g))
  {
    Point_3 p = get(vpm, vd);
    IO::internal::output_properties(os, &p, make_ply_point_writer (CGAL::Identity_property_map<Point_3>()));
    if(has_vcolor)
    {
      const CGAL::Color& c = get(vcm, vd);
      if(get_mode(os) == CGAL::IO::ASCII)
        os << c << std::endl;
      else
        os.write(reinterpret_cast<const char*>(&c), sizeof(c));
    }
  }

  std::vector<std::size_t> polygon;
  for(face_descriptor fd : faces(g))
  {
    polygon.clear();
    for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, g), g))
      polygon.push_back(get(vim, target(hd,g)));

    IO::internal::output_properties(os, &polygon,
                                    std::make_pair(CGAL::Identity_property_map<std::vector<std::size_t> >(),
                                                   PLY_property<std::vector<int> >("vertex_indices")));
    if(has_fcolor)
    {
      const CGAL::Color& c = get(fcm, fd);
      if(get_mode(os) == CGAL::IO::ASCII)
        os << c << std::endl;
      else
        os.write(reinterpret_cast<const char*>(&c), sizeof(c));
    }
  }

  return os.good();
}

template <typename Graph>
bool write_PLY(std::ostream& os, const Graph& g, const std::string& comments,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_PLY(os, g, comments, parameters::all_default());
}

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(std::ostream& os, const Graph& g, const CGAL_BGL_NP_CLASS& np,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_PLY(os, g, "", np);
}

template <typename Graph>
bool write_PLY(std::ostream& os, const Graph& g,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_PLY(os, g, "", parameters::all_default());
}

/*!
  \ingroup PkgBGLIOFct

  writes the graph in the output file `fname` in PLY format.

 \tparam Graph a model of `FaceListGraph`
 \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

 \param fname the name of the output file
 \param g the graph to be output
 \param comments a string included line by line in the header of the PLY stream (each line will be precedeed by "comment ")
 \param np optional \ref bgl_namedparameters "Named Parameters" described below

 \cgalNamedParamsBegin
    \cgalParamNBegin{use_binary_mode}
      \cgalParamDescription{indicates whether data should be written in binary (`true`) or in ASCII (`false`)}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd

   \cgalParamNBegin{vertex_index_map}
     \cgalParamDescription{a property map associating to each vertex of `graph` a unique index between `0` and `num_vertices(graph) - 1`}
     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                    as key type and `std::size_t` as value type}
     \cgalParamDefault{no vertex indices in the output}
   \cgalParamNEnd

    \cgalParamNBegin{vertex_color_map}
      \cgalParamDescription{a property map associating colors to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{no vertex color in the output}
    \cgalParamNEnd

    \cgalParamNBegin{face_color_map}
      \cgalParamDescription{a property map associating colors to the faces of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%face_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{no face color in the output}
    \cgalParamNEnd
 \cgalNamedParamsEnd

 \returns `true` if writing was successful.

 \see \ref IOStreamPLY
*/
template <typename Graph, typename NamedParameters>
bool write_PLY(const char* fname,
               const Graph& g,
               const std::string& comments,
               const NamedParameters& np
#ifndef DOXYGEN_RUNNING
              , typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
               )
{
  std::ofstream os(fname);
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
    CGAL::set_mode(os, CGAL::IO::BINARY);

  return write_PLY(os, g, comments, np);
}

template <typename Graph>
bool write_PLY(const char* fname, const Graph& g, const std::string comments,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_PLY(fname, g, comments, parameters::all_default());
}

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(const char* fname, const Graph& g, const CGAL_BGL_NP_CLASS& np,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_PLY(fname, g, "", np);
}

template <typename Graph>
bool write_PLY(const char* fname, const Graph& g,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_PLY(fname, g, "", parameters::all_default());
}

template <typename Graph, typename NamedParameters>
bool write_PLY(const std::string& fname, const Graph& g, const std::string& comments, const NamedParameters& np,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_PLY(fname.c_str(), g, comments, np);
}

template <typename Graph>
bool write_PLY(const std::string& fname, const Graph& g, const std::string comments,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_PLY(fname, g, comments, parameters::all_default());
}

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_PLY(const std::string& fname, const Graph& g, const CGAL_BGL_NP_CLASS& np,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_PLY(fname, g, "", np);
}

template <typename Graph>
bool write_PLY(const std::string& fname, const Graph& g,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_PLY(fname, g, "", parameters::all_default());
}

} // namespace CGAL

#endif // CGAL_BGL_IO_PLY_H
