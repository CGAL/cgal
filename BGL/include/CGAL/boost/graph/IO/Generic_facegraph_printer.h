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
#include <CGAL/property_map.h>

#include <boost/container/flat_map.hpp>

#include <string>
#include <vector>

namespace CGAL{
namespace IO {
namespace internal {

template <typename Stream, typename FaceGraph, typename FileWriter>
class Generic_facegraph_printer
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor                   vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type                  vertices_size_type;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor                     face_descriptor;

public:
  Generic_facegraph_printer(Stream& out) : m_out(out) { }
  Generic_facegraph_printer(Stream& out, FileWriter writer) : m_out(out), m_writer(writer) { }

  template <typename NamedParameters>
  bool operator()(const FaceGraph& g,
                  const NamedParameters& np)
  {
    typedef typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::const_type   VPM;
    typedef typename boost::property_traits<VPM>::reference                            Point_ref;

    typedef typename Polygon_mesh_processing::GetK<FaceGraph, NamedParameters>::Kernel Kernel;
    typedef typename Kernel::Vector_3                                                  Vector;
    typedef typename Kernel::Point_2                                                   Texture;
    typedef CGAL::Color                                                                Color;

    typedef typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::type         VPM;

    typedef typename internal_np::Lookup_named_param_def<
      internal_np::vertex_normal_map_t, NamedParameters,
      Constant_property_map<vertex_descriptor, Vector> >::type                         VNM;
    typedef typename internal_np::Lookup_named_param_def<
      internal_np::vertex_color_map_t, NamedParameters,
      Constant_property_map<vertex_descriptor, Color> >::type                          VCM;
    typedef typename internal_np::Lookup_named_param_def<
      internal_np::vertex_texture_map_t, NamedParameters,
      Constant_property_map<vertex_descriptor, Texture> >::type                        VTM;
    typedef typename internal_np::Lookup_named_param_def<
      internal_np::face_color_map_t, NamedParameters,
      Constant_property_map<face_descriptor, Color> >::type                            FCM;

    using parameters::choose_parameter;
    using parameters::is_default_parameter;
    using parameters::get_parameter;

    if(!m_out.good())
      return false;

    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, g));

    const bool has_vertex_normals = !(is_default_parameter(get_parameter(np, internal_np::vertex_normal_map)));
    const bool has_vertex_colors = !(is_default_parameter(get_parameter(np, internal_np::vertex_color_map)));
    const bool has_vertex_textures = !(is_default_parameter(get_parameter(np, internal_np::vertex_texture_map)));
    const bool has_face_colors = !(is_default_parameter(get_parameter(np, internal_np::face_color_map)));

    VNM vnm = choose_parameter(get_parameter(np, internal_np::vertex_normal_map), VNM());
    VCM vcm = choose_parameter(get_parameter(np, internal_np::vertex_color_map), VCM());
    VTM vtm = choose_parameter(get_parameter(np, internal_np::vertex_texture_map), VTM());
    FCM fcm = choose_parameter(get_parameter(np, internal_np::face_color_map), FCM());

    // @todo bench that against CGAL::Inverse_index and std::unordered_map
    boost::container::flat_map<vertex_descriptor, vertices_size_type> index_map;

    m_writer.write_header(m_out, num_vertices(g), num_halfedges(g), num_faces(g));

    vertices_size_type id = 0;
    for(const vertex_descriptor v : vertices(g))
    {
      const Point_ref p = get(vpm, v);
      m_writer.write_vertex(::CGAL::to_double(p.x()), ::CGAL::to_double(p.y()), ::CGAL::to_double(p.z()));

      if(has_vertex_normals)
      {
        const Vector& n = get(vnm, v);
        m_writer.write_vertex_normal(::CGAL::to_double(n.x()), ::CGAL::to_double(n.y()), ::CGAL::to_double(n.z()));
      }

      if(has_vertex_colors)
      {
        const CGAL::Color& vc = get(vcm, v);
        m_writer.write_vertex_color(vc.red(), vc.green(), vc.blue()); // @fixme correct?
      }

      if(has_vertex_textures)
      {
        const Texture& t = get(vtm, v);
        m_writer.write_vertex_texture(::CGAL::to_double(t.x()), ::CGAL::to_double(t.y()));
      }

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

      if(has_face_colors)
      {
        const CGAL::Color& fc = get(fcm, f);
        m_writer.write_face_normal(fc.red(), fc.green(), fc.blue());
      }

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
