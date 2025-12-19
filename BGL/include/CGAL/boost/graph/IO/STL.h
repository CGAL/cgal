// Copyright (c) 2015-2020  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//                 Mael Rouxel-Labbé

#ifndef CGAL_BGL_IO_STL_H
#define CGAL_BGL_IO_STL_H

#include <CGAL/boost/graph/IO/Generic_facegraph_builder.h>
#include <CGAL/IO/STL.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <fstream>
#include <iostream>
#include <string>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {
namespace internal {

// Use CRTP to gain access to the protected members without getters/setters.
template <typename Graph, typename Point>
class STL_builder
  : public Generic_facegraph_builder<Graph, Point, STL_builder<Graph, Point> >
{
  typedef STL_builder<Graph, Point>                                             Self;
  typedef Generic_facegraph_builder<Graph, Point, Self>                         Base;

  typedef typename Base::Point_container                                        Point_container;
  typedef typename Base::Face                                                   Face;
  typedef typename Base::Face_container                                         Face_container;

public:
  STL_builder(std::istream& is) : Base(is) { }

  template <typename NamedParameters>
  bool read(std::istream& is,
            Point_container& points,
            Face_container& faces,
            const NamedParameters& np)
  {
    return read_STL(is, points, faces, np);
  }
};

} // namespace internal

/*!
  \ingroup PkgBGLIoFuncsSTL

  \brief reads the graph `g` from the input stream, using the \ref IOStreamSTL.

  The data is expected to represent a 2-manifold (possibly with borders).

  \attention The graph `g` is not cleared, and the data from the stream are appended.

  \attention To read a binary file, the flag `std::ios::binary` must be set during the creation of the `ifstream`.

  \tparam Graph a model of `MutableFaceGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param is the input stream
  \param g the graph to be built from the input data
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd

    \cgalParamNBegin{verbose}
      \cgalParamDescription{whether extra information is printed when an incident occurs during reading}
      \cgalParamType{Boolean}
      \cgalParamDefault{`false`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if reading was successful and the resulting mesh is valid, `false` otherwise.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
*/
template <typename Graph, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_STL(std::istream& is,
              Graph& g,
              const CGAL_NP_CLASS& np)
{
  typedef typename CGAL::GetVertexPointMap<Graph, CGAL_NP_CLASS>::type      VPM;
  typedef typename boost::property_traits<VPM>::value_type                      Point;
  if(!is.good())
    return false;
  internal::STL_builder<Graph, Point> builder(is);
  return builder(g, np);
}

/*!
  \ingroup PkgBGLIoFuncsSTL

  \brief reads the graph `g` from the file `fname`, using the \ref IOStreamSTL.

  The data is expected to represent a 2-manifold (possibly with borders).
  If `use_binary_mode` is `true`, but the reading fails, ASCII reading will be automatically tested.

  \attention The graph `g` is not cleared, and the data from the file are appended.

  \tparam Graph a model of `MutableFaceGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the name of the input file
  \param g the graph to be built from the input data
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamNBegin{use_binary_mode}
      \cgalParamDescription{indicates whether data should be read in binary (`true`) or in \ascii (`false`)}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd

    \cgalParamNBegin{verbose}
      \cgalParamDescription{whether extra information is printed when an incident occurs during reading}
      \cgalParamType{Boolean}
      \cgalParamDefault{`false`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if reading was successful and the resulting mesh is valid, `false` otherwise.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
*/
template <typename Graph, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_STL(const std::string& fname,
              Graph& g, const
              CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  const bool binary = choose_parameter(get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ifstream is(fname, std::ios::binary);
    CGAL::IO::set_mode(is, CGAL::IO::BINARY);
    if(read_STL(is, g, np))
    {
      return true;
    }
    clear(g);
  }
  std::ifstream is(fname);
  CGAL::IO::set_mode(is, CGAL::IO::ASCII);

  typedef typename CGAL::GetVertexPointMap<Graph, CGAL_NP_CLASS>::type      VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(CGAL::vertex_point, g));
  bool v = choose_parameter(get_parameter(np, internal_np::verbose),
                            false);
  return read_STL(is, g, CGAL::parameters::use_binary_mode(false).vertex_point_map(vpm).verbose(v));
}

/// \cond SKIP_IN_MANUAL

template <typename Graph>
bool read_STL(std::istream& is, Graph& g) { return read_STL(is, g, parameters::default_values()); }

/// \endcond

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
  \ingroup PkgBGLIoFuncsSTL

  \brief writes the graph `g` in the output stream `os`, using the \ref IOStreamSTL.

  \attention To write to a binary file, the flag `std::ios::binary` must be set during the creation
             of the `ofstream`, and the \link PkgStreamSupportEnumRef `IO::Mode` \endlink of the stream
             must be set to `BINARY`.

  \tparam Graph a model of `FaceListGraph` and `HalfedgeListGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param os the output stream
  \param g the graph to be written
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

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
       \cgalParamType{int}
       \cgalParamDefault{the precision of the stream `os`}
       \cgalParamExtra{This parameter is only meaningful while using \ascii encoding.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre The graph must contain only triangle faces.

  \returns `true` if writing was successful, `false` otherwise.
*/
template <typename Graph, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_STL(std::ostream& os,
               const Graph& g,
               const CGAL_NP_CLASS& np = parameters::default_values())
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor                  halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor                      face_descriptor;

  typedef typename CGAL::GetVertexPointMap<Graph, CGAL_NP_CLASS>::const_type        VPM;
  typedef typename boost::property_traits<VPM>::reference                           Point_ref;

  typedef typename GetGeomTraits<Graph, CGAL_NP_CLASS>::type                        Kernel;
  typedef typename Kernel::Vector_3                                                 Vector;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  Kernel k = choose_parameter<Kernel>(get_parameter(np, internal_np::geom_traits));

  if(!os.good())
    return false;

  set_stream_precision_from_NP(os, np);

  if(get_mode(os) == BINARY)
  {
    os << "FileType: Binary                                                                ";
    const std::uint32_t N32 = static_cast<std::uint32_t>(faces(g).size());
    os.write(reinterpret_cast<const char *>(&N32), sizeof(N32));

    for(const face_descriptor f : faces(g))
    {
      const halfedge_descriptor h = halfedge(f, g);
      Point_ref p = get(vpm, target(h, g));
      Point_ref q = get(vpm, target(next(h, g), g));
      Point_ref r = get(vpm, source(h, g));

      const Vector n = internal::construct_normal_of_STL_face(p, q, r, k);

      const float coords[12] =
      {
        static_cast<float>(to_double(n.x())), static_cast<float>(to_double(n.y())), static_cast<float>(to_double(n.z())),
        static_cast<float>(to_double(p.x())), static_cast<float>(to_double(p.y())), static_cast<float>(to_double(p.z())),
        static_cast<float>(to_double(q.x())), static_cast<float>(to_double(q.y())), static_cast<float>(to_double(q.z())),
        static_cast<float>(to_double(r.x())), static_cast<float>(to_double(r.y())), static_cast<float>(to_double(r.z())) };

      for(int i=0; i<12; ++i)
        os.write(reinterpret_cast<const char *>(&coords[i]), sizeof(coords[i]));
      os << "  ";
    }
    os << std::flush;
  }
  else
  {
    os << "solid" << std::endl;
    for(const face_descriptor f : faces(g))
    {
      const halfedge_descriptor h = halfedge(f, g);
      Point_ref p = get(vpm, target(h, g));
      Point_ref q = get(vpm, target(next(h, g), g));
      Point_ref r = get(vpm, source(h, g));
      const Vector n = internal::construct_normal_of_STL_face(p, q, r, k);

      os << "facet normal " << n << "\nouter loop"<< "\n";
      os << "vertex " << p << "\n";
      os << "vertex " << q << "\n";
      os << "vertex " << r << "\n";
      os << "endloop\nendfacet" << "\n";
    }
    os << "endsolid" << std::endl;
  }

  return os.good();
}

/*!
  \ingroup PkgBGLIoFuncsSTL

  \brief writes the graph `g` into a file named `fname`, using the \ref IOStreamSTL.

  \tparam Graph a model of `FaceListGraph` and `HalfedgeListGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the name of the output stream
  \param g the graph to be written
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamNBegin{use_binary_mode}
      \cgalParamDescription{indicates whether data should be written in binary (`true`) or in \ascii (`false`)}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
       \cgalParamType{int}
       \cgalParamDefault{`6`}
       \cgalParamExtra{This parameter is only meaningful while using \ascii encoding.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre The graph must contain only triangle faces.

  \returns `true` if writing was successful, `false` otherwise.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
*/
template <typename Graph, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_STL(const std::string& fname, const Graph& g, const CGAL_NP_CLASS& np = parameters::default_values())
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ofstream os(fname, std::ios::binary);
    CGAL::IO::set_mode(os, CGAL::IO::BINARY);
    return write_STL(os, g, np);
  }
  else
  {
    std::ofstream os(fname);
    CGAL::IO::set_mode(os, CGAL::IO::ASCII);

    return write_STL(os, g, np);
  }
}

}} // namespace CGAL::IO

#endif // CGAL_BGL_IO_STL_H
