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
/// Read

namespace IO {
namespace internal {

// Use CRTP to gain access to the protected members without getters/setters.
template <typename FaceGraph, typename Point>
class STL_builder
  : public Generic_facegraph_builder<FaceGraph, Point, STL_builder<FaceGraph, Point> >
{
  typedef STL_builder<FaceGraph, Point>                                         Self;
  typedef Generic_facegraph_builder<FaceGraph, Point, Self>                     Base;

  typedef typename Base::Point_container                                        Point_container;
  typedef typename Base::Face                                                   Face;
  typedef typename Base::Face_container                                         Face_container;

public:
  STL_builder(std::istream& is_, bool verbose) : Base(is_, verbose) { }

  template <typename NamedParameters>
  bool read(std::istream& input,
            Point_container& points,
            Face_container& faces,
            const NamedParameters& np,
            bool verbose)
  {
    return read_STL(input, points, faces, np, verbose);
  }
};

} // namespace internal
} // namespace IO

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the STL format.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \pre The data must represent a 2-manifold

  \see \ref IOStreamSTL
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_STL(std::istream& in,
              FaceGraph& g,
              const CGAL_BGL_NP_CLASS& np,
              bool verbose = true)
{
  typedef typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  IO::internal::STL_builder<FaceGraph, Point> builder(in, verbose);
  return builder(g, np);
}

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from the file `fname` in the STL format.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \pre The data must represent a 2-manifold

  \see \ref IOStreamSTL
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_STL(const char* fname, FaceGraph& g, const CGAL_BGL_NP_CLASS& np,
              bool verbose = true)
{
  std::ifstream in(fname);
  return read_STL(in, g, np, verbose);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_STL(const std::string& fname, FaceGraph& g, const CGAL_BGL_NP_CLASS& np,
              bool verbose = true)
{
  return read_STL(fname.c_str(), g, np, verbose);
}

template <typename FaceGraph>
bool read_STL(std::istream& is, FaceGraph& g) { return read_STL(is, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_STL(const char* fname, FaceGraph& g) { return read_STL(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_STL(const std::string& fname, FaceGraph& g) { return read_STL(fname, g, parameters::all_default()); }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write


/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the stream `out` in the STL format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
  \cgalNamedParamsEnd
  \pre The graph must contain only triangle faces.

  \see \ref IOStreamSTL
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_STL(std::ostream& out,
               const FaceGraph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor                  halfedge_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor                      face_descriptor;

  typedef typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::const_type    VPM;
  typedef typename boost::property_traits<VPM>::reference                               Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                              Point;
  typedef typename Kernel_traits<Point>::Kernel::Vector_3                               Vector;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  if(!out.good())
    return false;

  if(get_mode(out) == IO::BINARY)
  {
    out << "FileType: Binary                                                                ";
    const boost::uint32_t N32 = static_cast<boost::uint32_t>(faces(g).size());
    out.write(reinterpret_cast<const char *>(&N32), sizeof(N32));

    for(const face_descriptor f : faces(g))
    {
      const halfedge_descriptor h = halfedge(f, g);
      Point_ref p = get(vpm, target(h, g));
      Point_ref q = get(vpm, target(next(h, g), g));
      Point_ref r = get(vpm, source(h, g));

      Vector n = collinear(p, q, r) ? Vector(1, 0, 0) : unit_normal(p, q, r);

      const float coords[12] =
      {
        static_cast<float>(to_double(n.x())), static_cast<float>(to_double(n.y())), static_cast<float>(to_double(n.z())),
        static_cast<float>(to_double(p.x())), static_cast<float>(to_double(p.y())), static_cast<float>(to_double(p.z())),
        static_cast<float>(to_double(q.x())), static_cast<float>(to_double(q.y())), static_cast<float>(to_double(q.z())),
        static_cast<float>(to_double(r.x())), static_cast<float>(to_double(r.y())), static_cast<float>(to_double(r.z())) };

      for(int i=0; i<12; ++i)
        out.write(reinterpret_cast<const char *>(&coords[i]), sizeof(coords[i]));
      out << "  ";
    }
  }
  else
  {
    out << "solid\n";
    for(const face_descriptor f : faces(g))
    {
      halfedge_descriptor h = halfedge(f, g);
      Point_ref p = get(vpm, target(h, g));
      Point_ref q = get(vpm, target(next(h, g), g));
      Point_ref r = get(vpm, source(h, g));
      Vector n = collinear(p, q, r) ? Vector(1, 0, 0) : unit_normal(p, q, r);

      out << "facet normal " << n << "\nouter loop\n";
      out << "vertex " << p << "\n";
      out << "vertex " << q << "\n";
      out << "vertex " << r << "\n";
      out << "endloop\nendfacet\n";
    }
    out << "endsolid\n";
  }

  return out.good();
}

/*!
\ingroup PkgBGLIOFct

 writes the graph `g` in the STL format into a file named `fname`.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
  \cgalNamedParamsEnd

 \sa Overloads of this function for specific models of the concept `FaceGraph`.

 \see \ref IOStreamSTL
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_STL(const char* fname, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream out(fname);
  return write_STL(out, g, np);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_STL(const std::string& fname, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  return write_STL(fname.c_str(), g, np);
}

template <typename FaceGraph>
bool write_STL(std::ostream& os, const FaceGraph& g) { return write_STL(os, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_STL(const char* fname, const FaceGraph& g) { return write_STL(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_STL(const std::string& fname, const FaceGraph& g) { return write_STL(fname, g, parameters::all_default()); }

} // namespace CGAL

#endif // CGAL_BGL_IO_STL_H
