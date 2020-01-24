// Copyright (c) 2019 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_BGL_IO_GENERIC_FACEGRAPH_BUILDER_H
#define CGAL_BGL_IO_GENERIC_FACEGRAPH_BUILDER_H

#include <CGAL/boost/graph/Euler_operations.h>

#include <iostream>
#include <string>
#include <vector>

namespace CGAL{
namespace IO {
namespace internal {

template <typename FaceGraph, typename Point, typename Derived>
class Generic_facegraph_builder
{
protected:
  typedef std::vector<Point>                                                    Point_container;
  typedef typename Point_container::size_type                                   size_type;
  typedef std::vector<std::size_t>                                              Face;
  typedef std::vector<Face>                                                     Face_container;

  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor            vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor              face_descriptor;

public:
  Generic_facegraph_builder(std::istream& in_) : m_in(in_) { }

  // must be implemented by the derived class
  bool read()
  {
    if(!m_in.good())
      return false;

    return static_cast<Derived*>(this)->read(m_in, m_points, m_faces);
  }

  template <typename NamedParameters>
  bool construct(FaceGraph& g,
                 const NamedParameters& np)
  {
    typedef typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::type  VPM;
    VPM vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                           get_property_map(CGAL::vertex_point, g));

    std::vector<vertex_descriptor> vertices(m_points.size());
    for(std::size_t id = 0, ps=m_points.size(); id<ps; ++id)
    {
      vertices[id] = add_vertex(g);
      put(vpm, vertices[id], m_points[id]);
    }

    for(size_type i=0, fs=m_faces.size(); i<fs; ++i)
    {
      std::vector<vertex_descriptor> face(m_faces[i].size());
      for(std::size_t j=0, fs=face.size(); j<fs; ++j)
        face[j] = vertices[m_faces[i][j]];

      face_descriptor f = CGAL::Euler::add_face(face, g);
      if(f == boost::graph_traits<FaceGraph>::null_face())
        return false;
    }

    return is_valid_polygon_mesh(g);
  }

  template <typename NamedParameters>
  bool operator()(FaceGraph& g, const NamedParameters& np)
  {
    if(!read())
      return false;

    return construct(g, np);
  }

  std::string name, color;

protected:
  std::istream& m_in;

  Point_container m_points;
  Face_container m_faces;
};

} // end internal
} // end IO
} // end CGAL

#endif // CGAL_BGL_IO_GENERIC_FACEGRAPH_BUILDER_H
