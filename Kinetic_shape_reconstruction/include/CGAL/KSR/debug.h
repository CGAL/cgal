// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_DEBUG_H
#define CGAL_KSR_DEBUG_H

#if defined(WIN32) || defined(_WIN32)
#define _NL_ "\r\n"
#else
#define _NL_ "\n"
#endif

// STL includes.
#include <fstream>

// CGAL includes.
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Random.h>

// Internal includes.
#include <CGAL/KSR/utils.h>

namespace CGAL {
namespace KSR_3 {

std::tuple<unsigned char, unsigned char, unsigned char>
get_idx_color (KSR::size_t idx) {

  CGAL::Random rand (idx);
  return std::make_tuple ((unsigned char)(rand.get_int(32, 192)),
                          (unsigned char)(rand.get_int(32, 192)),
                          (unsigned char)(rand.get_int(32, 192)));
}

template <typename DS>
void dump_intersection_edges (const DS& data, const std::string& tag = std::string()) {

  std::string filename = (tag != std::string() ? tag + "_" : "") + "intersection_edges.polylines.txt";
  std::ofstream out (filename);
  out.precision(18);

  for (const typename DS::IEdge iedge : data.iedges())
    out << "2 " << data.segment_3 (iedge) << std::endl;
}

template <typename DS>
void dump_segmented_edges (const DS& data, const std::string& tag = std::string()) {

  std::vector<std::ofstream*> out;
  for (KSR::size_t i = 0; i < data.nb_intersection_lines(); ++ i) {
    std::string filename = (tag != std::string() ? tag + "_" : "") + "intersection_line_" + std::to_string(i) + ".polylines.txt";
    out.push_back (new std::ofstream (filename));
    out.back()->precision(18);
  }

  for (const typename DS::IEdge iedge : data.iedges()) {
    CGAL_assertion (data.line_idx(iedge) != KSR::no_element());
    *(out[data.line_idx(iedge)]) << "2 " << data.segment_3 (iedge) << std::endl;
  }

  for (std::ofstream* o : out)
    delete o;
}

template <typename DS>
void dump_constrained_edges (const DS& data, const std::string& tag = std::string()) {

  std::string filename = (tag != std::string() ? tag + "_" : "") + "constrained_edges.polylines.txt";
  std::ofstream out (filename);
  out.precision(18);

  for (KSR::size_t i = 0; i < data.number_of_support_planes(); ++ i) {
    for (const typename DS::PEdge pedge : data.pedges(i))
      if (data.has_iedge(pedge))
        out << "2 " << data.segment_3 (pedge) << std::endl;
  }
}

template <typename DS>
void dump_polygons (const DS& data, const std::string& tag = std::string()) {

  typedef CGAL::Surface_mesh<typename DS::Kernel::Point_3> Mesh;
  typedef typename Mesh::template Property_map<typename Mesh::Face_index, unsigned char> Uchar_map;

  Mesh mesh;
  Uchar_map red = mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("red", 0).first;
  Uchar_map green = mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("green", 0).first;
  Uchar_map blue = mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("blue", 0).first;

#ifdef CGAL_KSR_DEBUG

  Mesh dbg_mesh;
  Uchar_map dbg_red = dbg_mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("red", 0).first;
  Uchar_map dbg_green = dbg_mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("green", 0).first;
  Uchar_map dbg_blue = dbg_mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("blue", 0).first;

#endif

  Mesh bbox_mesh;
  Uchar_map bbox_red = bbox_mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("red", 0).first;
  Uchar_map bbox_green = bbox_mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("green", 0).first;
  Uchar_map bbox_blue = bbox_mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("blue", 0).first;

  KSR::vector<typename Mesh::Vertex_index> vertices;
  KSR::vector<typename Mesh::Vertex_index> map_vertices;

  for (KSR::size_t i = 0; i < data.number_of_support_planes(); ++ i) {
    if (data.is_bbox_support_plane(i)) {

      map_vertices.clear();
      for (typename DS::PVertex pvertex : data.pvertices(i)) {
        if (map_vertices.size() <= pvertex.second)
          map_vertices.resize (pvertex.second + 1);
        map_vertices[pvertex.second] = bbox_mesh.add_vertex (data.point_3(pvertex));
      }

      for (typename DS::PFace pface : data.pfaces(i)) {
        vertices.clear();
        for(typename DS::PVertex pvertex : data.pvertices_of_pface(pface))
          vertices.push_back (map_vertices[pvertex.second]);

        typename Mesh::Face_index face = bbox_mesh.add_face (vertices);
        std::tie (bbox_red[face], bbox_green[face], bbox_blue[face])
          = get_idx_color ((i+1) * (pface.second+1));
      }

    } else {

      map_vertices.clear();
      for (typename DS::PVertex pvertex : data.pvertices(i)) {
        if (map_vertices.size() <= pvertex.second)
          map_vertices.resize (pvertex.second + 1);
        map_vertices[pvertex.second] = mesh.add_vertex (data.point_3 (pvertex));
      }

      for (typename DS::PFace pface : data.pfaces(i)) {
        vertices.clear();
        for(typename DS::PVertex pvertex : data.pvertices_of_pface(pface))
          vertices.push_back (map_vertices[pvertex.second]);
        typename Mesh::Face_index face = mesh.add_face (vertices);
        std::tie (red[face], green[face], blue[face])
          = get_idx_color (i * (pface.second+1));
      }

#ifdef CGAL_KSR_DEBUG

      map_vertices.clear();
      for (typename DS::PVertex pvertex : data.dbg_pvertices(i)) {
        if (map_vertices.size() <= pvertex.second)
          map_vertices.resize (pvertex.second + 1);
        map_vertices[pvertex.second] = dbg_mesh.add_vertex (data.dbg_point_3 (pvertex));
      }

      for (typename DS::PFace pface : data.dbg_pfaces(i)) {
        vertices.clear();
        for(typename DS::PVertex pvertex : data.dbg_pvertices_of_pface(pface))
          vertices.push_back (map_vertices[pvertex.second]);
        typename Mesh::Face_index face = dbg_mesh.add_face (vertices);
        std::tie (dbg_red[face], dbg_green[face], dbg_blue[face])
          = get_idx_color (i * (pface.second+1));
      }

#endif

    }
  }

  std::string filename = (tag != std::string() ? tag + "_" : "") + "polygons.ply";
  std::ofstream out (filename);
  // CGAL::set_binary_mode (out);
  CGAL::write_ply(out, mesh);

#ifdef CGAL_KSR_DEBUG

  std::string dbg_filename = (tag != std::string() ? tag + "_" : "") + "dbg_polygons.ply";
  std::ofstream dbg_out (dbg_filename);
  // CGAL::set_binary_mode (dbg_out);
  CGAL::write_ply(dbg_out, dbg_mesh);

#endif

#if 0

  std::string bbox_filename = (tag != std::string() ? tag + "_" : "") + "bbox_polygons.ply";
  std::ofstream bbox_out (bbox_filename);
  // CGAL::set_binary_mode (bbox_out);
  CGAL::write_ply(bbox_out, bbox_mesh);

#endif

}

template <typename DS>
void dump_polygon_borders (const DS& data, const std::string& tag = std::string()) {

  std::string filename = (tag != std::string() ? tag + "_" : "") + "polygon_borders.polylines.txt";
  std::ofstream out (filename);

  for (KSR::size_t i = 6; i < data.number_of_support_planes(); ++ i)
    for (const typename DS::PEdge pedge : data.pedges(i))
      out << "2 " << data.segment_3 (pedge) << std::endl;

  // {
  //   std::string filename = (tag != std::string() ? tag + "_" : "") + "polygon_borders_perturbated.polylines.txt";
  //   std::ofstream out (filename);

  //   CGAL::Random r;
  //   for (KSR::size_t i = 6; i < data.number_of_support_planes(); ++ i)
  //     for (const typename DS::PEdge pedge : data.pedges(i))
  //     {
  //       typename DS::Kernel::Point_3 s = data.segment_3 (pedge).source ();
  //       s = s + typename DS::Kernel::Vector_3 (r.get_double(-0.01, 0.01),r.get_double(-0.01, 0.01),r.get_double(-0.01, 0.01));
  //       typename DS::Kernel::Point_3 t = data.segment_3 (pedge).target ();
  //       CGAL::Random rt (t.x() * t.y() * t.z());
  //       t = t + typename DS::Kernel::Vector_3 (r.get_double(-0.01, 0.01),r.get_double(-0.01, 0.01),r.get_double(-0.01, 0.01));
  //       out << "2 " <<  s << " " << t << std::endl;
  //     }
  // }
}

template <typename DS, typename Event>
void dump_event (const DS& data, const Event& ev, const std::string& tag = std::string()) {

  if (ev.is_pvertex_to_pvertex()) {

    std::string vfilename = (tag != std::string() ? tag + "_" : "") + "event_pvertex.xyz";
    std::ofstream vout (vfilename);
    vout.precision(18);
    vout << data.point_3 (ev.pvertex()) << std::endl;

    std::string ofilename = (tag != std::string() ? tag + "_" : "") + "event_pother.xyz";
    std::ofstream oout (ofilename);
    oout.precision(18);
    oout << data.point_3 (ev.pother()) << std::endl;

  } else if (ev.is_pvertex_to_iedge()) {

    std::string lfilename = (tag != std::string() ? tag + "_" : "") + "event_iedge.polylines.txt";
    std::ofstream lout (lfilename);
    lout.precision(18);
    lout << "2 " << data.segment_3 (ev.iedge()) << std::endl;

    std::string vfilename = (tag != std::string() ? tag + "_" : "") + "event_pvertex.xyz";
    std::ofstream vout (vfilename);
    vout.precision(18);
    vout << data.point_3 (ev.pvertex());

  } else if (ev.is_pvertex_to_ivertex()) {

    std::string vfilename = (tag != std::string() ? tag + "_" : "") + "event_pvertex.xyz";
    std::ofstream vout (vfilename);
    vout.precision(18);
    vout << data.point_3 (ev.pvertex()) << std::endl;

    std::string ofilename = (tag != std::string() ? tag + "_" : "") + "event_ivertex.xyz";
    std::ofstream oout (ofilename);
    oout.precision(18);
    oout << data.point_3 (ev.ivertex()) << std::endl;
  }
}

template <typename DS>
void dump (const DS& data, const std::string& tag = std::string()) {

  dump_intersection_edges (data, tag);
  // dump_constrained_edges (data, tag);
  // dump_polygon_borders (data, tag);
  dump_polygons (data, tag);
}

template<typename GeomTraits>
class Saver {

public:
  using Traits    = GeomTraits;
  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Point_3   = typename Traits::Point_3;
  using Segment_2 = typename Traits::Segment_2;
  using Segment_3 = typename Traits::Segment_3;

  using Color        = CGAL::Color;
  using Surface_mesh = CGAL::Surface_mesh<Point_3>;
  using Random       = CGAL::Random;

  Saver() :
  m_path_prefix("/Users/monet/Documents/gf/kinetic/logs/"),
  grey(Color(125, 125, 125)),
  red(Color(125, 0, 0))
  { }

  void initialize(std::stringstream& stream) const {
    stream.precision(20);
  }

  void export_points_2(
    const KSR::vector<Point_2>& points,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    for (const auto& point : points)
      stream << point << " 0 " << std::endl;
    save(stream, file_name + ".xyz");
  }

  void export_points_3(
    const KSR::vector<Point_3>& points,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    for (const auto& point : points)
      stream << point << std::endl;
    save(stream, file_name + ".xyz");
  }

  void export_segments_2(
    const KSR::vector<Segment_2>& segments,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    for (const auto& segment : segments)
      stream << "2 " << segment.source() << " 0 " << segment.target() << " 0 " << std::endl;
    save(stream, file_name + ".polylines");
  }

  void export_segments_3(
    const KSR::vector<Segment_3>& segments,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    for (const auto& segment : segments)
      stream << "2 " << segment.source() << " " << segment.target() << std::endl;
    save(stream, file_name + ".polylines");
  }

  void export_polygon_soup_3(
    const std::vector< std::vector<Point_3> >& polygons,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    std::size_t num_vertices = 0;
    for (const auto& polygon : polygons)
      num_vertices += polygon.size();
    std::size_t num_faces = polygons.size();
    add_ply_header_mesh(stream, num_vertices, num_faces);

    for (const auto& polygon : polygons)
      for (const auto& p : polygon)
        stream << p << std::endl;

    std::size_t i = 0, polygon_id = 0;
    for (const auto& polygon : polygons) {
      stream << polygon.size() << " ";
      for (std::size_t j = 0; j < polygon.size(); ++j)
        stream << i++ << " ";
      stream << get_idx_color(polygon_id) << std::endl;
      ++polygon_id;
    }
    save(stream, file_name + ".ply");
  }

  void export_polygon_soup_3(
    const KSR::vector< KSR::vector<Point_3> >& polygons,
    const KSR::vector<Color>& colors,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    KSR::size_t num_vertices = 0;
    for (const auto& polygon : polygons)
      num_vertices += polygon.size();
    KSR::size_t num_faces = polygons.size();
    add_ply_header_mesh(stream, num_vertices, num_faces);

    for (const auto& polygon : polygons)
      for (const auto& p : polygon)
        stream << p << std::endl;

    KSR::size_t i = 0, polygon_id = 0;
    for (KSR::size_t k = 0; k < polygons.size(); ++k) {
      stream << polygons[k].size() << " ";
      for (KSR::size_t j = 0; j < polygons[k].size(); ++j)
        stream << i++ << " ";
      stream << colors[k] << std::endl;
      ++polygon_id;
    }
    save(stream, file_name + ".ply");
  }

  void export_bounding_box_3(
    const KSR::array<Point_3, 8>& bounding_box,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    Surface_mesh bbox;
    CGAL_assertion(bounding_box.size() == 8);
    CGAL::make_hexahedron(
      bounding_box[0], bounding_box[1], bounding_box[2], bounding_box[3],
      bounding_box[4], bounding_box[5], bounding_box[6], bounding_box[7], bbox);
    stream << bbox;
    save(stream, file_name + ".off");
  }

  void export_mesh_2(
    const KSR::vector<Point_2>& vertices,
    const KSR::vector< KSR::vector<KSR::size_t> >& faces,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    add_ply_header_mesh(stream, vertices.size(), faces.size());
    for (const auto& vertex : vertices)
      stream << vertex << " 0 " << std::endl;

    for (const auto& face : faces) {
      stream << face.size();
      for (const KSR::size_t findex : face){
        stream << " " << findex;
      }
      stream << " " << grey << std::endl;
    }
    save(stream, file_name + ".ply");
  }

  void export_mesh_2(
    const KSR::vector<Point_2>& vertices,
    const KSR::vector< KSR::vector<KSR::size_t> >& faces,
    const KSR::vector<Color>& colors,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    add_ply_header_mesh(stream, vertices.size(), faces.size());
    for (const auto& vertex : vertices)
      stream << vertex << " 0 " << std::endl;

    for (KSR::size_t k = 0; k < faces.size(); ++k) {
      stream << faces[k].size();
      for (const KSR::size_t findex : faces[k]){
        stream << " " << findex;
      }
      stream << " " << colors[k] << std::endl;
    }
    save(stream, file_name + ".ply");
  }

  const Color get_idx_color(const KSR::size_t idx) const {
    Random rand(idx);
    const unsigned char r = rand.get_int(32, 192);
    const unsigned char g = rand.get_int(32, 192);
    const unsigned char b = rand.get_int(32, 192);
    return Color(r, g, b);
  }

private:
  const std::string m_path_prefix;
  const Color grey, red;

  void save(
    const std::stringstream& stream,
    const std::string file_name) const {

    const std::string file_path = m_path_prefix + file_name;
    std::ofstream file(file_path.c_str(), std::ios_base::out);

    if (!file)
      std::cerr << std::endl <<
        "ERROR: Error saving file " << file_path
      << "!" << std::endl << std::endl;

    file << stream.str();
    file.close();
  }

  void add_ply_header_points(
    std::stringstream& stream,
    const KSR::size_t size) const {

    stream <<
    "ply" 				         +  std::string(_NL_) + ""               			 <<
    "format ascii 1.0"     +  std::string(_NL_) + ""     			           <<
    "element vertex "      << size        << "" + std::string(_NL_) + "" <<
    "property double x"    +  std::string(_NL_) + ""    			           <<
    "property double y"    +  std::string(_NL_) + ""    			           <<
    "property double z"    +  std::string(_NL_) + "" 				             <<
    "property uchar red"   +  std::string(_NL_) + "" 				             <<
    "property uchar green" +  std::string(_NL_) + "" 				             <<
    "property uchar blue"  +  std::string(_NL_) + "" 				             <<
    "property uchar alpha" +  std::string(_NL_) + "" 				             <<
    "end_header"           +  std::string(_NL_) + "";
  }

  void add_ply_header_normals(
    std::stringstream& stream,
    const KSR::size_t size) const {

    stream <<
    "ply" 				         +  std::string(_NL_) + ""               			 <<
    "format ascii 1.0"     +  std::string(_NL_) + ""     			           <<
    "element vertex "      << size        << "" + std::string(_NL_) + "" <<
    "property double x"    +  std::string(_NL_) + ""    			           <<
    "property double y"    +  std::string(_NL_) + ""    			           <<
    "property double z"    +  std::string(_NL_) + "" 				             <<
    "property double nx"   +  std::string(_NL_) + ""    			           <<
    "property double ny"   +  std::string(_NL_) + ""    			           <<
    "property double nz"   +  std::string(_NL_) + "" 				             <<
    "end_header"           +  std::string(_NL_) + "";
  }

  void add_ply_header_mesh(
    std::stringstream& stream,
    const KSR::size_t num_vertices,
    const KSR::size_t num_faces) const {

    stream <<
    "ply" 				         +  std::string(_NL_) + ""               			<<
    "format ascii 1.0"     +  std::string(_NL_) + ""     			          <<
    "element vertex "      << num_vertices     << "" + std::string(_NL_) + "" <<
    "property double x"    +  std::string(_NL_) + ""    			          <<
    "property double y"    +  std::string(_NL_) + ""    			          <<
    "property double z"    +  std::string(_NL_) + "" 				            <<
    "element face "        << num_faces        << "" + std::string(_NL_) + "" <<
    "property list uchar int vertex_indices"         + std::string(_NL_) + "" <<
    "property uchar red"   +  std::string(_NL_) + "" 				            <<
    "property uchar green" +  std::string(_NL_) + "" 				            <<
    "property uchar blue"  +  std::string(_NL_) + "" 				            <<
    "property uchar alpha" +  std::string(_NL_) + "" 				             <<
    "end_header"           +  std::string(_NL_) + "";
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_DEBUG_H
