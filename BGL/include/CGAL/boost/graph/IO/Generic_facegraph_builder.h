// Copyright (c) 2019 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno
//                 Mael Rouxel-Labbé

#ifndef CGAL_BGL_IO_GENERIC_FACEGRAPH_BUILDER_H
#define CGAL_BGL_IO_GENERIC_FACEGRAPH_BUILDER_H

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>
#include <string>
#include <vector>

namespace CGAL{
namespace IO {
namespace internal {

template <typename Graph, typename Point, typename Derived>
class Generic_facegraph_builder
{
protected:
  typedef std::vector<Point>                                                           Point_container;
  typedef typename Point_container::size_type                                          size_type;
  typedef std::vector<std::size_t>                                                     Face;
  typedef std::vector<Face>                                                            Face_container;

  typedef typename boost::graph_traits<Graph>::vertex_descriptor                       vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor                         face_descriptor;

public:
  Generic_facegraph_builder(std::istream& in_) : m_is(in_) { }

  template <typename NamedParameters = parameters::Default_named_parameters>
  bool operator()(Graph& g, const NamedParameters& np = parameters::default_values())
  {
    typedef typename GetK<Graph, NamedParameters>::Kernel                              Kernel;
    typedef typename Kernel::Vector_3                                                  Vector;
    typedef typename Kernel::Point_2                                                   Texture;
    typedef CGAL::IO::Color                                                            Color;

    typedef typename CGAL::GetVertexPointMap<Graph, NamedParameters>::type             VPM;

    // usually will be true, but might not be the case if using custom type points
//    static_assert(std::is_same<typename Kernel::Point_3,
//                                        typename boost::property_traits<VPM>::value_type>::value);
//    static_assert(std::is_same<typename Kernel::Point_3,
//                                        typename boost::range_value<Point_container>::type>::value);

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

    typedef typename boost::property_traits<VNM>::value_type                           Vertex_normal;
    typedef typename boost::property_traits<VCM>::value_type                           Vertex_color;
    typedef typename boost::property_traits<VTM>::value_type                           Vertex_texture;
    typedef typename boost::property_traits<FCM>::value_type                           Face_color;

    using parameters::choose_parameter;
    using parameters::is_default_parameter;
    using parameters::get_parameter;

    const bool is_vnm_requested = !(is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value);
    const bool is_vcm_requested = !(is_default_parameter<NamedParameters, internal_np::vertex_color_map_t>::value);
    const bool is_vtm_requested = !(is_default_parameter<NamedParameters, internal_np::vertex_texture_map_t>::value);
    const bool is_fcm_requested = !(is_default_parameter<NamedParameters, internal_np::face_color_map_t>::value);

    std::vector<Vertex_normal> vertex_normals;
    std::vector<Vertex_color> vertex_colors;
    std::vector<Vertex_texture> vertex_textures;
    std::vector<Face_color> face_colors;

    const bool verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);
    const bool binary = choose_parameter(get_parameter(np, internal_np::use_binary_mode), true);

    bool ok =
        static_cast<Derived*>(this)->read(m_is, m_points, m_faces,
                                          parameters::vertex_normal_output_iterator(std::back_inserter(vertex_normals))
                                                     .vertex_color_output_iterator(std::back_inserter(vertex_colors))
                                                     .vertex_texture_output_iterator(std::back_inserter(vertex_textures))
                                                     .face_color_output_iterator(std::back_inserter(face_colors))
                                                     .verbose(verbose)
                                                     .use_binary_mode(binary));
    if(!ok)
      return false;

    // Construct the graph
    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point), get_property_map(CGAL::vertex_point, g));
    VNM vnm = choose_parameter(get_parameter(np, internal_np::vertex_normal_map), VNM());
    VCM vcm = choose_parameter(get_parameter(np, internal_np::vertex_color_map), VCM());
    VTM vtm = choose_parameter(get_parameter(np, internal_np::vertex_texture_map), VTM());
    FCM fcm = choose_parameter(get_parameter(np, internal_np::face_color_map), FCM());

    const bool has_vertex_normals = (is_vnm_requested && !(vertex_normals.empty()));
    const bool has_vertex_colors = (is_vcm_requested && !(vertex_colors.empty()));
    const bool has_vertex_textures = (is_vtm_requested && !(vertex_textures.empty()));
    const bool has_face_colors = (is_fcm_requested && !(face_colors.empty()));

    if(has_vertex_normals && vertex_normals.size() != m_points.size())
      return false;
    if(has_vertex_colors && vertex_colors.size() != m_points.size())
      return false;
    if(has_vertex_textures && vertex_textures.size() != m_points.size())
      return false;
    if(has_face_colors && face_colors.size() != m_faces.size())
      return false;

    std::vector<vertex_descriptor> vertices(m_points.size());

    for(std::size_t id=0, ps=m_points.size(); id<ps; ++id)
    {
      vertices[id] = add_vertex(g);
      put(vpm, vertices[id], m_points[id]);

      // extra properties
      if(has_vertex_normals)
        put(vnm, vertices[id], vertex_normals[id]);
      if(has_vertex_colors)
        put(vcm, vertices[id], vertex_colors[id]);
      if(has_vertex_textures)
        put(vtm, vertices[id], vertex_textures[id]);
    }

    for(size_type i=0, fs=m_faces.size(); i<fs; ++i)
    {
      std::vector<vertex_descriptor> face(m_faces[i].size());
      for(std::size_t j=0, fis=face.size(); j<fis; ++j)
        face[j] = vertices[m_faces[i][j]];

      face_descriptor f = CGAL::Euler::add_face(face, g);
      if(f == boost::graph_traits<Graph>::null_face())
        return false;

      if(has_face_colors)
        put(fcm, f, face_colors[i]);
    }

    return is_valid(g);
  }

protected:
  std::istream& m_is;

  Point_container m_points;
  Face_container m_faces;
};

} // namespace internal
} // namespace IO
} // namespace CGAL

#endif // CGAL_BGL_IO_GENERIC_FACEGRAPH_BUILDER_H
