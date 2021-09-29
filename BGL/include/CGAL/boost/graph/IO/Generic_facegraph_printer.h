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

// Unfortunately, we don't know the value type of the normal/texture property maps
template <typename VNM>
struct Normal_writer
{
  Normal_writer(const VNM vnm) : vnm(vnm) { }

  template <typename Writer, typename VD>
  void operator()(Writer& writer, const VD v) const
  {
    const typename boost::property_traits<VNM>::reference n = get(vnm, v);
    writer.write_vertex_normal(to_double(n.x()), to_double(n.y()), to_double(n.z()));
  }

private:
  const VNM vnm;
};

template <>
struct Normal_writer<internal_np::Param_not_found>
{
  Normal_writer(const internal_np::Param_not_found&) { }

  template <typename Writer, typename VD>
  void operator()(Writer&, const VD) const { }
};

template <typename VTM>
struct Texture_writer
{
  Texture_writer(const VTM vtm) : vtm(vtm) { }

  template <typename Writer, typename VD>
  void operator()(Writer& writer, const VD v) const
  {
    const typename boost::property_traits<VTM>::reference t = get(vtm, v);
    writer.write_vertex_texture(to_double(t.x()), to_double(t.y()));
  }

private:
  const VTM vtm;
};

template <>
struct Texture_writer<internal_np::Param_not_found>
{
  Texture_writer(const internal_np::Param_not_found&) { }

  template <typename Writer, typename VD>
  void operator()(Writer&, const VD) const { }
};

template <typename Stream, typename Graph, typename FileWriter>
class Generic_facegraph_printer
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor                   vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::vertices_size_type                  vertices_size_type;
  typedef typename boost::graph_traits<Graph>::face_descriptor                     face_descriptor;

public:
  Generic_facegraph_printer(Stream& os) : m_os(os) { }
  Generic_facegraph_printer(Stream& os, FileWriter writer) : m_os(os), m_writer(writer) { }

  template <typename NamedParameters>
  bool operator()(const Graph& g,
                  const NamedParameters& np)
  {
    typedef typename GetVertexPointMap<Graph, NamedParameters>::const_type         VPM;
    typedef typename boost::property_traits<VPM>::reference                        Point_ref;

    typedef CGAL::IO::Color                                                        Color;

    typedef typename internal_np::Lookup_named_param_def<
      internal_np::vertex_color_map_t, NamedParameters,
      Constant_property_map<vertex_descriptor, Color> >::type                      VCM;
    typedef typename internal_np::Lookup_named_param_def<
      internal_np::face_color_map_t, NamedParameters,
      Constant_property_map<face_descriptor, Color> >::type                        FCM;

    // No default because value_type is unknown, but the pmap is only used if provided via NP
    typedef typename internal_np::Get_param<
      typename NamedParameters::base, internal_np::vertex_normal_map_t>::type      VNM;
    typedef typename internal_np::Get_param<
      typename NamedParameters::base, internal_np::vertex_texture_map_t>::type     VTM;

    using parameters::choose_parameter;
    using parameters::is_default_parameter;
    using parameters::get_parameter;

    if(!m_os.good())
      return false;

    set_stream_precision_from_NP(m_os, np);

    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, g));

    const bool has_vertex_normals = !(is_default_parameter(get_parameter(np, internal_np::vertex_normal_map)));
    const bool has_vertex_colors = !(is_default_parameter(get_parameter(np, internal_np::vertex_color_map)));
    const bool has_vertex_textures = !(is_default_parameter(get_parameter(np, internal_np::vertex_texture_map)));
    const bool has_face_colors = !(is_default_parameter(get_parameter(np, internal_np::face_color_map)));

    VNM vnm = get_parameter(np, internal_np::vertex_normal_map);
    VTM vtm = get_parameter(np, internal_np::vertex_texture_map);
    VCM vcm = choose_parameter<VCM>(get_parameter(np, internal_np::vertex_color_map));
    FCM fcm = choose_parameter<FCM>(get_parameter(np, internal_np::face_color_map));

    Normal_writer<VNM> nw(vnm);
    Texture_writer<VTM> tw(vtm);

    // @todo bench that against CGAL::Inverse_index and std::unordered_map
    boost::container::flat_map<vertex_descriptor, vertices_size_type> index_map;
    m_writer.write_header(m_os, vertices(g).size(), halfedges(g).size(), faces(g).size(),
                          has_face_colors || has_vertex_colors,
                          has_vertex_normals                  ,
                          has_vertex_textures                 );

    vertices_size_type id = 0;
    for(const vertex_descriptor v : vertices(g))
    {
      const Point_ref p = get(vpm, v);
      m_writer.write_vertex(to_double(p.x()), to_double(p.y()), to_double(p.z()));

      if(has_vertex_normals)
        nw(m_writer, v);

      if(has_vertex_colors)
      {
        const CGAL::IO::Color& vc = get(vcm, v);
        m_writer.write_vertex_color(vc.red(), vc.green(), vc.blue()); // @fixme correct?
      }

      if(has_vertex_textures)
        tw(m_writer, v);

      index_map[v] = id++;
    }

    m_writer.write_facet_header();
    for(const face_descriptor f : faces(g))
    {
      CGAL::Halfedge_around_face_circulator<Graph> hc(halfedge(f, g), g);
      CGAL::Halfedge_around_face_circulator<Graph> hc_end = hc;

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
        const CGAL::IO::Color& fc = get(fcm, f);
        m_writer.write_face_color(fc.red(), fc.green(), fc.blue());
      }

      m_writer.write_facet_end();
    }
    m_writer.write_footer();

    return m_os.good();
  }

  bool operator()(const Graph& g) { return operator()(g, parameters::all_default()); }

protected:
  Stream& m_os;
  FileWriter m_writer;
};

} // end internal
} // end IO
} // end CGAL

#endif // CGAL_BGL_IO_GENERIC_FACEGRAPH_PRINTER_H
