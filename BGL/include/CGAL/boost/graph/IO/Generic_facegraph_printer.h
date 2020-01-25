// Copyright (c) 2019 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_BGL_IO_GENERIC_FACEGRAPH_PRINTER_H
#define CGAL_BGL_IO_GENERIC_FACEGRAPH_PRINTER_H

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/container/flat_map.hpp>

#include <string>
#include <vector>

namespace CGAL{
namespace IO {
namespace internal {

template <typename Stream, typename FaceGraph, typename FileWriter>
class Generic_facegraph_printer
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor                  vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type                 vertices_size_type;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor                    face_descriptor;

public:
  Generic_facegraph_printer(Stream& out) : m_out(out) { }
  Generic_facegraph_printer(Stream& out, FileWriter writer) : m_out(out), m_writer(writer) { }

  template <typename NamedParameters>
  bool operator()(const FaceGraph& g,
                  const NamedParameters& np)
  {
    typedef typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::const_type  VPM;
    typedef typename boost::property_traits<VPM>::reference                           Point_ref;

    VPM vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                           get_const_property_map(CGAL::vertex_point, g));

    if(!m_out.good())
      return false;

    // @todo bench that against CGAL::Inverse_index and std::unordered_map
    boost::container::flat_map<vertex_descriptor, vertices_size_type> index_map;

    m_writer.write_header(m_out, num_vertices(g), num_halfedges(g), num_faces(g));

    vertices_size_type id = 0;
    for(const vertex_descriptor v : vertices(g))
    {
      const Point_ref p = get(vpm, v);
      m_writer.write_vertex(::CGAL::to_double(p.x()),
                            ::CGAL::to_double(p.y()),
                            ::CGAL::to_double(p.z()));
      index_map[v] = id++;
    }

    m_writer.write_facet_header();
    for(const face_descriptor f : faces(g))
    {
      CGAL::Halfedge_around_face_circulator<FaceGraph> hc(halfedge(f, g), g);
      CGAL::Halfedge_around_face_circulator<FaceGraph> hc_end = hc;

      const std::size_t n = circulator_size(hc);
      CGAL_assertion(n >= 3);

      m_writer.write_facet_begin(n);
      do
      {
        m_writer.write_facet_vertex_index(index_map[target(*hc, g)]);
        ++hc;
      }
      while(hc != hc_end);

      m_writer.write_facet_end();
    }
    m_writer.write_footer();

    return m_out.good();
  }

  bool operator()(const FaceGraph& g) { return operator()(g, parameters::all_default()); }

protected:
  Stream& m_out;
  FileWriter m_writer;
};

} // end internal
} // end IO
} // end CGAL

#endif // CGAL_BGL_IO_GENERIC_FACEGRAPH_PRINTER_H
