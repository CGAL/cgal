// Copyright (c) 2023 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Florent Lafarge, Dmitry Anisimov, Simon Giraudot

#ifndef CGAL_KSP_DEBUG_H
#define CGAL_KSP_DEBUG_H

#include <CGAL/license/Kinetic_space_partition.h>

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
#include <CGAL/Surface_mesh/IO.h>
#include <CGAL/Random.h>
#include <CGAL/IO/Color.h>
#include <CGAL/boost/graph/generators.h>

// Internal includes.
#include <CGAL/KSP/utils.h>

namespace CGAL {
namespace KSP_3 {
namespace internal {

const std::tuple<unsigned char, unsigned char, unsigned char>
get_idx_color(std::size_t idx) {

  CGAL::Random rand(static_cast<unsigned int>(idx));
  return std::make_tuple(
    static_cast<unsigned char>(rand.get_int(32, 192)),
    static_cast<unsigned char>(rand.get_int(32, 192)),
    static_cast<unsigned char>(rand.get_int(32, 192)));
}

template<typename DS>
void dump_intersection_edges(const DS& data, const std::string tag = std::string()) {

  const std::string filename = (tag != std::string() ? tag + "-" : "") + "intersection-edges.polylines.txt";
  std::ofstream out(filename);
  out.precision(20);

  for (const auto iedge : data.iedges()) {
    out << "2 " << data.segment_3(iedge) << std::endl;
  }
  out.close();
}

template<typename DS>
void dump_segmented_edges(const DS& data, const std::string tag = std::string()) {

  std::vector<std::ofstream*> out;
  for (std::size_t i = 0; i < data.nb_intersection_lines(); ++i) {
    const std::string filename = (tag != std::string() ? tag + "-" : "") + "intersection-line-" + std::to_string(i) + ".polylines.txt";
    out.push_back(new std::ofstream(filename));
    out.back()->precision(20);
  }

  for (const auto iedge : data.iedges()) {
    CGAL_assertion(data.line_idx(iedge) != std::size_t(-1));
    *(out[data.line_idx(iedge)]) << "2 " << data.segment_3(iedge) << std::endl;
  }

  for (std::ofstream* o : out) {
    delete o;
  }
}

template<typename DS>
void dump_constrained_edges(const DS& data, const std::string tag = std::string()) {

  const std::string filename = (tag != std::string() ? tag + "-" : "") + "constrained-edges.polylines.txt";
  std::ofstream out(filename);
  out.precision(20);

  for (std::size_t i = 0; i < data.number_of_support_planes(); ++i) {
    for (const auto pedge : data.pedges(i)) {
      if (data.has_iedge(pedge)) {
        out << "2 " << data.segment_3(pedge) << std::endl;
      }
    }
  }
  out.close();
}

template<typename DS>
void dump_2d_surface_mesh(
  const DS& data,
  const std::size_t support_plane_idx,
  const std::string tag = std::string()) {

  using Point_3      = typename DS::Kernel::Point_3;
  using Mesh         = CGAL::Surface_mesh<Point_3>;
  using Face_index   = typename Mesh::Face_index;
  using Vertex_index = typename Mesh::Vertex_index;
  using Int_map    = typename Mesh::template Property_map<Face_index, int>;

  Mesh mesh;
  Int_map red   = mesh.template add_property_map<Face_index, int>("red", 0).first;
  Int_map green = mesh.template add_property_map<Face_index, int>("green", 0).first;
  Int_map blue  = mesh.template add_property_map<Face_index, int>("blue", 0).first;

  std::vector<Vertex_index> vertices;
  std::vector<Vertex_index> map_vertices;

  map_vertices.clear();

  for (const auto pvertex : data.pvertices(support_plane_idx)) {
    if (map_vertices.size() <= pvertex.second) {
      map_vertices.resize(pvertex.second + 1);
    }
    //pts.push_back(data.point_3(pvertex));
    map_vertices[pvertex.second] = mesh.add_vertex(data.point_3(pvertex));
  }

  for (const auto pface : data.pfaces(support_plane_idx)) {
    vertices.clear();
    for (const auto pvertex : data.pvertices_of_pface(pface)) {
      vertices.push_back(map_vertices[pvertex.second]);
    }

    CGAL_assertion(vertices.size() >= 3);
    const auto face = mesh.add_face(vertices);
    if (face == Mesh::null_face()) {
      std::cout << "could not dump mesh " << tag << std::endl;
      return;
    }
    std::tie(red[face], green[face], blue[face]) =
      get_idx_color((support_plane_idx + 1) * (pface.second + 1));
  }

  const std::string filename = (tag != std::string() ? tag + "-" : "") + "polygons.ply";
  std::ofstream out(filename);
  out.precision(20);
  CGAL::IO::write_PLY(out, mesh);
  out.close();
}

template<typename DS>
void dump_polygons(const DS& data, const std::string tag = std::string()) {

  using Point_3      = typename DS::Point_3;
  using Mesh         = CGAL::Surface_mesh<Point_3>;
  using Face_index   = typename Mesh::Face_index;
  using Vertex_index = typename Mesh::Vertex_index;
  using Uchar_map    = typename Mesh::template Property_map<Face_index, unsigned char>;

  Mesh mesh;
  Uchar_map red   = mesh.template add_property_map<Face_index, unsigned char>("red", 0).first;
  Uchar_map green = mesh.template add_property_map<Face_index, unsigned char>("green", 0).first;
  Uchar_map blue  = mesh.template add_property_map<Face_index, unsigned char>("blue", 0).first;

  Mesh bbox_mesh;
  Uchar_map bbox_red   = bbox_mesh.template add_property_map<Face_index, unsigned char>("red", 0).first;
  Uchar_map bbox_green = bbox_mesh.template add_property_map<Face_index, unsigned char>("green", 0).first;
  Uchar_map bbox_blue  = bbox_mesh.template add_property_map<Face_index, unsigned char>("blue", 0).first;

  std::vector<Vertex_index> vertices;
  std::vector<Vertex_index> map_vertices;

  for (std::size_t i = 0; i < data.number_of_support_planes(); ++i) {
    if (data.is_bbox_support_plane(i)) {

      map_vertices.clear();
      for (const auto pvertex : data.pvertices(i)) {
        if (map_vertices.size() <= pvertex.second) {
          map_vertices.resize(pvertex.second + 1);
        }
        map_vertices[pvertex.second] = bbox_mesh.add_vertex(data.point_3(pvertex));
      }

      for (const auto pface : data.pfaces(i)) {
        vertices.clear();
        for (const auto pvertex : data.pvertices_of_pface(pface)) {
          vertices.push_back(map_vertices[pvertex.second]);
        }

        CGAL_assertion(vertices.size() >= 3);
        const auto face = bbox_mesh.add_face(vertices);
        CGAL_assertion(face != Mesh::null_face());
        std::tie(bbox_red[face], bbox_green[face], bbox_blue[face]) =
          get_idx_color((i + 1) * (pface.second + 1));
      }

    } else {

      map_vertices.clear();
      for (const auto pvertex : data.pvertices(i)) {
        if (map_vertices.size() <= pvertex.second) {
          map_vertices.resize(pvertex.second + 1);
        }
        map_vertices[pvertex.second] = mesh.add_vertex(data.point_3(pvertex));
      }

      for (const auto pface : data.pfaces(i)) {
        vertices.clear();
        for (const auto pvertex : data.pvertices_of_pface(pface)) {
          vertices.push_back(map_vertices[pvertex.second]);
        }

        CGAL_assertion(vertices.size() >= 3);
        const auto face = mesh.add_face(vertices);
        CGAL_assertion(face != Mesh::null_face());
        std::tie(red[face], green[face], blue[face]) =
          get_idx_color(i * (pface.second + 1));
      }
    }
  }

  const std::string filename = (tag != std::string() ? tag + "-" : "") + "polygons.ply";
  std::ofstream out(filename);
  out.precision(20);
  CGAL::IO::write_PLY(out, mesh);
  out.close();

#if false

  const std::string bbox_filename = (tag != std::string() ? tag + "-" : "") + "bbox-polygons.ply";
  std::ofstream bbox_out(bbox_filename);
  bbox_out.precision(20);
  CGAL::write_PLY(bbox_out, bbox_mesh);
  bbox_out.close();

#endif
}

template<typename DS>
void dump_polygon_borders(const DS& data, const std::string tag = std::string()) {

  const std::string filename = (tag != std::string() ? tag + "-" : "") + "polygon-borders.polylines.txt";
  std::ofstream out(filename);
  out.precision(20);
  for (std::size_t i = 6; i < data.number_of_support_planes(); ++i) {
    for (const auto pedge : data.pedges(i)) {
      out << "2 " << data.segment_3(pedge) << std::endl;
    }
  }
  out.close();
}

template<typename DS, typename Event>
void dump_event(const DS& data, const Event& ev, const std::string tag = std::string()) {

  if (ev.is_pvertex_to_pvertex()) {

    const std::string vfilename = (tag != std::string() ? tag + "-" : "") + "event-pvertex.xyz";
    std::ofstream vout(vfilename);
    vout.precision(20);
    vout << data.point_3(ev.pvertex()) << std::endl;
    vout.close();

    const std::string ofilename = (tag != std::string() ? tag + "-" : "") + "event-pother.xyz";
    std::ofstream oout(ofilename);
    oout.precision(20);
    oout << data.point_3(ev.pother()) << std::endl;
    oout.close();

  } else if (ev.is_pvertex_to_iedge()) {

    const std::string lfilename = (tag != std::string() ? tag + "-" : "") + "event-iedge.polylines.txt";
    std::ofstream lout(lfilename);
    lout.precision(20);
    lout << "2 " << data.segment_3(ev.iedge()) << std::endl;
    lout.close();

    const std::string vfilename = (tag != std::string() ? tag + "-" : "") + "event-pvertex.xyz";
    std::ofstream vout(vfilename);
    vout.precision(20);
    vout << data.point_3(ev.pvertex());
    vout.close();

  } else if (ev.is_pvertex_to_ivertex()) {

    const std::string vfilename = (tag != std::string() ? tag + "-" : "") + "event-pvertex.xyz";
    std::ofstream vout(vfilename);
    vout.precision(20);
    vout << data.point_3(ev.pvertex()) << std::endl;
    vout.close();

    const std::string ofilename = (tag != std::string() ? tag + "-" : "") + "event-ivertex.xyz";
    std::ofstream oout(ofilename);
    oout.precision(20);
    oout << data.point_3(ev.ivertex()) << std::endl;
    oout.close();
  }
}

template<typename DS>
void dump(const DS& data, const std::string tag = std::string()) {

  dump_polygons(data, tag);
  dump_intersection_edges(data, tag);
}

template<typename GeomTraits>
class Saver {

public:
  using Traits    = GeomTraits;
  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Point_3   = typename Traits::Point_3;
  using Vector_3  = typename Traits::Vector_3;
  using Segment_2 = typename Traits::Segment_2;
  using Segment_3 = typename Traits::Segment_3;

  using Color        = CGAL::IO::Color;
  using Surface_mesh = CGAL::Surface_mesh<Point_3>;
  using Random       = CGAL::Random;

  Saver() :
  m_path_prefix(""),
  grey(Color(125, 125, 125)),
  red(Color(125, 0, 0))
  { }

  void initialize(std::stringstream& stream) const {
    stream.precision(20);
  }

  void export_points_2(
    const std::vector<Point_2>& points,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    for (const auto& point : points)
      stream << point << " 0 " << std::endl;
    save(stream, file_name + ".xyz");
  }

  void export_points_2(
    const std::vector< std::vector<Point_2> >& regions,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    std::size_t num_points = 0;
    for (const auto& region : regions)
      num_points += region.size();
    add_ply_header_points(stream, num_points);

    for (std::size_t i = 0; i < regions.size(); ++i) {
      const auto color = get_idx_color(i);
      for (const auto& point : regions[i])
        stream << point << " 0 " << color << std::endl;
    }
    save(stream, file_name + ".ply");
  }

  void export_points_3(
    const std::vector<Point_3>& points,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    for (const auto& point : points)
      stream << point << std::endl;
    save(stream, file_name + ".xyz");
  }

  void export_regions_3(
    const std::vector<Point_3>& points,
    const std::vector<Vector_3>& normals,
    const std::vector<int>& region,
    const std::string file_name) const {

    if (points.size() != normals.size()) {
      std::cout << "export_regions_3: number of points and normals are not equal" << std::endl;
      return;
    }

    if (points.size() != region.size()) {
      std::cout << "export_regions_3: number of points and region indices are not equal" << std::endl;
      return;
    }

    std::stringstream stream;
    initialize(stream);

    add_ply_header_regions(stream, points.size());

    for (std::size_t i = 0; i < points.size(); ++i) {
      stream << points[i] << " " << normals[i] << " " << region[i] << std::endl;
    }
    save(stream, file_name);
  }

  void export_points_3(
    const std::vector<Point_3>& points,
    const std::vector<Vector_3>& normals,
    const std::string file_name) const {

    if (points.size() != normals.size()) {
      std::cout << "export_regions_3: number of points and normals are not equal" << std::endl;
      return;
    }

    std::stringstream stream;
    initialize(stream);

    add_ply_header_normals(stream, points.size());

    for (std::size_t i = 0; i < points.size(); ++i) {
      stream << points[i] << " " << normals[i] << " " << std::endl;
    }
    save(stream, file_name);
  }

  void export_points_3(
    const std::vector<Point_3>& points,
    const std::vector<Vector_3>& normals,
    const std::vector<Color>& colors,
    const std::string file_name) const {

    if (points.size() != normals.size() && points.size() != colors.size()) {
      std::cout << "export_regions_3: number of points and normals or colors are not equal" << std::endl;
      return;
    }

    std::stringstream stream;
    initialize(stream);

    add_ply_header_normals_colors(stream, points.size());

    for (std::size_t i = 0; i < points.size(); ++i) {
      stream << points[i] << " " << normals[i] << " " << colors[i] << std::endl;
    }
    save(stream, file_name);
  }

  void export_segments_2(
    const std::vector<Segment_2>& segments,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    for (const auto& segment : segments)
      stream << "2 " << segment.source() << " 0 " << segment.target() << " 0 " << std::endl;
    save(stream, file_name + ".polylines.txt");
  }

  void export_segments_3(
    const std::vector<Segment_3>& segments,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    for (const auto& segment : segments)
      stream << "2 " << segment.source() << " " << segment.target() << std::endl;
    save(stream, file_name + ".polylines.txt");
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
    const std::vector<std::vector< std::vector<Point_3> > >& polygons,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    std::size_t num_vertices = 0;
    std::size_t num_faces = 0;
    for (const auto& region : polygons) {
      num_faces += region.size();
      for (const auto& polygon : region)
        num_vertices += polygon.size();
    }
    add_ply_header_mesh(stream, num_vertices, num_faces);

    for (const auto& region : polygons)
      for(const auto& polygon : region)
        for (const auto& p : polygon)
         stream << p << std::endl;

    std::size_t i = 0, region_id = 0;
    for (const auto& region : polygons) {
      for (const auto& polygon : region) {
        stream << polygon.size() << " ";
        for (std::size_t j = 0; j < polygon.size(); ++j)
          stream << i++ << " ";
        stream << get_idx_color(region_id) << std::endl;
      }
      ++region_id;
    }
    save(stream, file_name + ".ply");
  }

  void export_indexed_triangles_3(
    const std::vector<Point_3>& vertices,
    const std::vector<std::size_t>& tris,
    const std::string file_name) const {
    assert((tris.size() % 3) == 0);

    std::stringstream stream;
    initialize(stream);

    add_ply_header_mesh_no_color(stream, vertices.size(), tris.size() / 3);

    for (const auto& v : vertices)
      stream << v << std::endl;

    for (std::size_t i = 0; i < (tris.size() - 2); i += 3)
      stream << "3 " << tris[i] << " " << tris[i + 1] << " " << tris[i + 2] << std::endl;

    save(stream, file_name + ".ply");
  }

  void export_indexed_polygons_3(
    const std::vector<Point_3>& vertices,
    const std::vector<std::vector<std::size_t> >& polys,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    add_ply_header_mesh_no_color(stream, vertices.size(), polys.size());

    for (const auto& v : vertices)
      stream << v << std::endl;

    for (const auto& poly : polys) {
      stream << poly.size();
      for (std::size_t i = 0; i < poly.size(); i++)
        stream << " " << poly[i];
      stream << std::endl;
    }

    save(stream, file_name + ".ply");
  }

  void export_polygon_soup_3(
    const std::vector< std::vector<Point_3> >& polygons,
    const std::vector<Color>& colors,
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

    std::size_t i = 0;
    for (std::size_t k = 0; k < polygons.size(); ++k) {
      stream << polygons[k].size() << " ";
      for (std::size_t j = 0; j < polygons[k].size(); ++j)
        stream << i++ << " ";
      stream << colors[k] << std::endl;
    }
    save(stream, file_name + ".ply");
  }

  void export_bounding_box_3(
    const std::array<Point_3, 8>& bounding_box,
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
    const std::vector<Point_2>& vertices,
    const std::vector< std::vector<std::size_t> >& faces,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    add_ply_header_mesh(stream, vertices.size(), faces.size());
    for (const auto& vertex : vertices)
      stream << vertex << " 0 " << std::endl;

    for (const auto& face : faces) {
      stream << face.size();
      for (const std::size_t findex : face){
        stream << " " << findex;
      }
      stream << " " << grey << std::endl;
    }
    save(stream, file_name + ".ply");
  }

  void export_mesh_2(
    const std::vector<Point_2>& vertices,
    const std::vector< std::vector<std::size_t> >& faces,
    const std::vector<Color>& colors,
    const std::string file_name) const {

    std::stringstream stream;
    initialize(stream);

    add_ply_header_mesh(stream, vertices.size(), faces.size());
    for (const auto& vertex : vertices)
      stream << vertex << " 0 " << std::endl;

    for (std::size_t k = 0; k < faces.size(); ++k) {
      stream << faces[k].size();
      for (const std::size_t findex : faces[k]){
        stream << " " << findex;
      }
      stream << " " << colors[k] << std::endl;
    }
    save(stream, file_name + ".ply");
  }

  const Color get_idx_color(const std::size_t idx) const {
    Random rand(static_cast<unsigned int>(idx));
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
        "ERROR: WHILE SAVING FILE " << file_path
      << "!" << std::endl << std::endl;

    file << stream.str();
    file.close();
  }

  void add_ply_header_points(
    std::stringstream& stream,
    const std::size_t size) const {

    stream <<
    "ply"                  +  std::string(_NL_) + ""                      <<
    "format ascii 1.0"     +  std::string(_NL_) + ""                      <<
    "element vertex "      << size        << "" + std::string(_NL_) + "" <<
    "property double x"    +  std::string(_NL_) + ""                     <<
    "property double y"    +  std::string(_NL_) + ""                     <<
    "property double z"    +  std::string(_NL_) + ""                      <<
    "property uchar red"   +  std::string(_NL_) + ""                      <<
    "property uchar green" +  std::string(_NL_) + ""                      <<
    "property uchar blue"  +  std::string(_NL_) + ""                      <<
    "property uchar alpha" +  std::string(_NL_) + ""                      <<
    "end_header"           +  std::string(_NL_) + "";
  }

  void add_ply_header_normals(
    std::stringstream& stream,
    const std::size_t size) const {

    stream <<
      "ply" + std::string(_NL_) + "" <<
      "format ascii 1.0" + std::string(_NL_) + "" <<
      "element vertex " << size << "" + std::string(_NL_) + "" <<
      "property double x" + std::string(_NL_) + "" <<
      "property double y" + std::string(_NL_) + "" <<
      "property double z" + std::string(_NL_) + "" <<
      "property double nx" + std::string(_NL_) + "" <<
      "property double ny" + std::string(_NL_) + "" <<
      "property double nz" + std::string(_NL_) + "" <<
      "end_header" + std::string(_NL_) + "";
  }

  void add_ply_header_normals_colors(
    std::stringstream& stream,
    const std::size_t size) const {

    stream <<
      "ply" + std::string(_NL_) + "" <<
      "format ascii 1.0" + std::string(_NL_) + "" <<
      "element vertex " << size << "" + std::string(_NL_) + "" <<
      "property double x" + std::string(_NL_) + "" <<
      "property double y" + std::string(_NL_) + "" <<
      "property double z" + std::string(_NL_) + "" <<
      "property double nx" + std::string(_NL_) + "" <<
      "property double ny" + std::string(_NL_) + "" <<
      "property double nz" + std::string(_NL_) + "" <<
      "property uchar red" + std::string(_NL_) + "" <<
      "property uchar green" + std::string(_NL_) + "" <<
      "property uchar blue" + std::string(_NL_) + "" <<
      "property uchar alpha" + std::string(_NL_) + "" <<
      "end_header" + std::string(_NL_) + "";
  }

  void add_ply_header_regions(
    std::stringstream& stream,
    const std::size_t size) const {

    stream <<
      "ply" << std::endl <<
      "format ascii 1.0" << std::endl <<
      "element vertex " << size << std::endl <<
      "property double x" << std::endl <<
      "property double y" << std::endl <<
      "property double z" << std::endl <<
      "property double nx" << std::endl <<
      "property double ny" << std::endl <<
      "property double nz" << std::endl <<
      "property int region" << std::endl <<
      "end_header" << std::endl;
  }

  void add_ply_header_mesh(
    std::stringstream& stream,
    const std::size_t num_vertices,
    const std::size_t num_faces) const {

    stream <<
      "ply" + std::string(_NL_) + "" <<
      "format ascii 1.0" + std::string(_NL_) + "" <<
      "element vertex " << num_vertices << "" + std::string(_NL_) + "" <<
      "property double x" + std::string(_NL_) + "" <<
      "property double y" + std::string(_NL_) + "" <<
      "property double z" + std::string(_NL_) + "" <<
      "element face " << num_faces << "" + std::string(_NL_) + "" <<
      "property list uchar int vertex_indices" + std::string(_NL_) + "" <<
      "property uchar red" + std::string(_NL_) + "" <<
      "property uchar green" + std::string(_NL_) + "" <<
      "property uchar blue" + std::string(_NL_) + "" <<
      "property uchar alpha" + std::string(_NL_) + "" <<
      "end_header" + std::string(_NL_) + "";
  }

  void add_ply_header_mesh_no_color(
    std::stringstream& stream,
    const std::size_t num_vertices,
    const std::size_t num_faces) const {

    stream <<
      "ply" + std::string(_NL_) + "" <<
      "format ascii 1.0" + std::string(_NL_) + "" <<
      "element vertex " << num_vertices << "" + std::string(_NL_) + "" <<
      "property double x" + std::string(_NL_) + "" <<
      "property double y" + std::string(_NL_) + "" <<
      "property double z" + std::string(_NL_) + "" <<
      "element face " << num_faces << "" + std::string(_NL_) + "" <<
      "property list ushort int vertex_indices" + std::string(_NL_) + "" <<
      "end_header" + std::string(_NL_) + "";
  }
};

template<typename DS, typename PFace>
void dump_volume(
  const DS& data,
  const std::vector<PFace>& pfaces,
  const std::string file_name,
  const bool use_colors = true,
  const std::size_t volume_index = 0) {

  using Point_3 = typename DS::Kernel::Point_3;
  std::vector<Point_3> polygon;
  std::vector< std::vector<Point_3> > polygons;
  std::vector<CGAL::IO::Color> colors;

  colors.reserve(pfaces.size());
  polygons.reserve(pfaces.size());

  Saver<typename DS::Kernel> saver;
  for (const auto& pface : pfaces) {
    const auto pvertices = data.pvertices_of_pface(pface);
    const auto color = saver.get_idx_color(volume_index);
    polygon.clear();
    for (const auto pvertex : pvertices) {
      polygon.push_back(data.point_3(pvertex));
    }
    if (use_colors) {
      colors.push_back(color);
    }
    else {
      colors.push_back(saver.get_idx_color(0));
    }
    CGAL_assertion(polygon.size() >= 3);
    polygons.push_back(polygon);
  }
  CGAL_assertion(colors.size() == pfaces.size());
  CGAL_assertion(polygons.size() == pfaces.size());

  saver.export_polygon_soup_3(polygons, colors, file_name);
}

template<typename DS, typename PFace, typename FT>
void dump_visi(
  const DS& data,
  const std::vector<PFace>& pfaces,
  const std::string &file_name,
  const FT color) {

  using Point_3 = typename DS::Kernel::Point_3;
  std::vector<Point_3> polygon;
  std::vector< std::vector<Point_3> > polygons;
  std::vector<CGAL::IO::Color> colors;

  colors.reserve(pfaces.size());
  polygons.reserve(pfaces.size());

  const CGAL::IO::Color low(255, 255, 255);
  const CGAL::IO::Color high(0, 0, 255);

  Saver<typename DS::Kernel> saver;
  for (const auto& pface : pfaces) {
    const auto pvertices = data.pvertices_of_pface(pface);
    polygon.clear();
    for (const auto pvertex : pvertices) {
      polygon.push_back(data.point_3(pvertex));
    }

    colors.push_back(CGAL::IO::Color(static_cast<unsigned char>((1 - color) * low[0] + color * high[0]), static_cast<unsigned char>((1 - color) * low[1] + color * high[1]), static_cast<unsigned char>((1 - color) * low[2] + color * high[2]), ((color > 0.5) ? 150 : 25)));

    CGAL_assertion(polygon.size() >= 3);
    polygons.push_back(polygon);
  }
  CGAL_assertion(colors.size() == pfaces.size());
  CGAL_assertion(polygons.size() == pfaces.size());

  saver.export_polygon_soup_3(polygons, colors, file_name);
}

template<typename DS>
void dump_volumes(const DS& data, const std::string tag = std::string()) {

  using Point_3 = typename DS::Kernel::Point_3;
  std::vector<Point_3> polygon;
  std::vector< std::vector<Point_3> > polygons;
  std::vector<CGAL::IO::Color> colors;

  Saver<typename DS::Kernel> saver;
  for (std::size_t i = 0; i < data.volumes().size(); ++i) {
    const auto& volume = data.volumes()[i];
    const auto color = saver.get_idx_color(i);

    colors.clear();
    polygons.clear();
    for (const auto& pface : volume.pfaces) {
      polygon.clear();
      for (const auto pvertex : data.pvertices_of_pface(pface)) {
        polygon.push_back(data.point_3(pvertex));
      }
      CGAL_assertion(polygon.size() >= 3);
      polygons.push_back(polygon);
      colors.push_back(color);
    }

    const std::string file_name =
      (tag != std::string() ? tag + "-" : "") + std::to_string(i);
    saver.export_polygon_soup_3(polygons, colors, file_name);
  }
}

template<typename K>
void dump_polygon(const std::vector<typename K::Point_3>& pts, const std::string& filename) {
  Saver<K> saver;
  std::vector<std::vector<typename K::Point_3> > pts2;
  pts2.push_back(pts);

  saver.export_polygon_soup_3(pts2, filename);
}

void dump_polygon(const std::vector<CGAL::Epick::Point_3>& pts, const std::string& filename) {
  Saver<CGAL::Epick> saver;
  std::vector<std::vector<CGAL::Epick::Point_3> > pts2;
  pts2.push_back(pts);

  saver.export_polygon_soup_3(pts2, filename);
}

void dump_polygona(const std::vector<CGAL::Epick::Point_3>& pts, const std::string& filename) {
  Saver<CGAL::Epick> saver;
  std::vector<std::vector<CGAL::Epick::Point_3> > pts2;
  pts2.push_back(pts);

  saver.export_polygon_soup_3(pts2, filename);
}

void dump_polygons(const std::vector<std::vector<CGAL::Epick::Point_3> >& pts, const std::string& filename) {
  Saver<CGAL::Epick> saver;

  saver.export_polygon_soup_3(pts, filename);
}
void dump_polygons(const std::vector<std::vector<std::vector<CGAL::Epick::Point_3> > >& pts, const std::string& filename) {
  Saver<CGAL::Epick> saver;

  saver.export_polygon_soup_3(pts, filename);
}

void dump_indexed_triangles(const std::vector<CGAL::Epick::Point_3>& pts, const std::vector<std::size_t>& tris, const std::string& filename) {
  Saver<CGAL::Epick> saver;

  saver.export_indexed_triangles_3(pts, tris, filename);
}

void dump_indexed_polygons(const std::vector<CGAL::Epick::Point_3>& pts, const std::vector<std::vector<std::size_t> >& polys, const std::string& filename) {
  Saver<CGAL::Epick> saver;

  saver.export_indexed_polygons_3(pts, polys, filename);
}

  void dump_polygons(const std::vector<std::vector<CGAL::Epick::Point_3> >& pts, const std::vector<CGAL::IO::Color>& colors, const std::string& filename) {
  Saver<CGAL::Epick> saver;

  saver.export_polygon_soup_3(pts, colors, filename);
}

template<typename DS, typename Polygon_2>
void dump_polygon(
  const DS& data,
  const std::size_t sp_idx,
  const Polygon_2& input,
  const std::string name) {

  using Kernel = typename DS::Kernel;
  using Point_3 = typename Kernel::Point_3;
  std::vector<Point_3> polygon;
  std::vector< std::vector<Point_3> > polygons;
  for (const auto& point_2 : input) {
    polygon.push_back(data.to_3d(sp_idx, point_2));
  }
  polygons.push_back(polygon);
  Saver<Kernel> saver;
  saver.export_polygon_soup_3(polygons, name);
}

template<typename DS, typename Polygon_2>
void dump_polygons(
  const DS& data,
  const Polygon_2& input,//    std::map< std::size_t, std::pair<Polygon_2, Indices> > polygons;
  const std::string name) {

  using Kernel = typename DS::Kernel;
  using Point_3 = typename Kernel::Point_3;
  std::vector<Point_3> polygon;
  std::vector< std::vector<Point_3> > polygons;
  for (const auto& pair : input) {
    for (const auto &point2d : pair.second.first)
      polygon.push_back(data.to_3d(pair.first, point2d));
    polygons.push_back(polygon);
    polygon.clear();
  }
  Saver<Kernel> saver;
  saver.export_polygon_soup_3(polygons, name);
}

void dump_points(const std::vector<CGAL::Epick::Point_3>& pts, const std::vector<CGAL::Epick::Vector_3>& normals, const std::vector<CGAL::IO::Color>& colors, const std::string& filename) {
  Saver<CGAL::Epick> saver;
  saver.export_points_3(pts, normals, colors, filename);
}

template<typename DS>
void dump_ifaces(const DS& data, const std::string tag = std::string()) {
  // write all polygons into a separate ply with support plane index and iface index
  for (std::size_t sp_idx = data.number_of_support_planes(); sp_idx++;) {
    for (typename DS::IFace f : data.support_plane(sp_idx).ifaces()) {
      Saver<typename DS::Kernel> saver;
      std::vector<std::vector<typename DS::Kernel::Point_3> > pts(1);
      for (auto v : data.igraph().face(f).vertices)
        pts.back().push_back(data.igraph().point_3(v));
      saver.export_polygon_soup_3(pts, tag + "-" + std::to_string(sp_idx) + "-" + std::to_string(f));
    }
  }
}
/*

template<typename DS, typename PTS>
void dump_polygon(const DS& data, PTS pts, std::string tag = std::string()) {
  // write all polygons into a separate ply with support plane index and iface index
  Saver<typename DS::Kernel> saver;
  saver.export_polygon_soup_3(pts, tag);
}*/

template<typename DS, typename PFace>
void dump_pface(
  const DS& data,
  const PFace& pface,
  const std::string name) {

  using Kernel = typename DS::Kernel;
  using Point_3 = typename Kernel::Point_3;
  std::vector<Point_3> polygon;
  std::vector< std::vector<Point_3> > polygons;
  for (const auto pvertex : data.pvertices_of_pface(pface)) {
    polygon.push_back(data.point_3(pvertex));
  }
  CGAL_assertion(polygon.size() >= 3);
  polygons.push_back(polygon);
  Saver<Kernel> saver;
  saver.export_polygon_soup_3(polygons, name);
}

template<typename DS, typename PFace>
void dump_pfaces(
  const DS& data,
  const std::vector<PFace>& pfaces,
  const std::string name) {

  using Kernel = typename DS::Kernel;
  using Point_3 = typename Kernel::Point_3;
  std::vector< std::vector<Point_3> > polygons;
  for (auto pface : pfaces) {
    std::vector<Point_3> polygon;
    for (const auto pvertex : data.pvertices_of_pface(pface)) {
      polygon.push_back(data.point_3(pvertex));
    }
    CGAL_assertion(polygon.size() >= 3);
    polygons.push_back(polygon);
  }
  Saver<Kernel> saver;
  saver.export_polygon_soup_3(polygons, name);
}

template<typename DS, typename PEdge>
void dump_pedge(
  const DS& data,
  const PEdge& pedge,
  const std::string name) {

  using Kernel = typename DS::Kernel;
  using Segment_3 = typename Kernel::Segment_3;
  const std::vector<Segment_3> segments = { data.segment_3(pedge) };
  Saver<Kernel> saver;
  saver.export_segments_3(segments, name);
}

template<typename DS, typename PFace, typename PEdge>
void dump_info(
  const DS& data,
  const PFace& pface,
  const PEdge& pedge,
  const std::vector<PFace>& nfaces,
  const std::string &suffix) {

  std::cout << "DEBUG: number of found nfaces: " << nfaces.size() << std::endl;
  dump_pface(data, pface, "face-curr-" + suffix);
  dump_pedge(data, pedge, "face-edge-" + suffix);
  for (std::size_t i = 0; i < nfaces.size(); ++i) {
    dump_pface(data, nfaces[i], "nface-" + std::to_string(i) + "-" + suffix);
  }
}

template<typename Point_3>
void dump_frame(
  const std::vector<Point_3>& points,
  const std::string name) {

  using Kernel = typename Kernel_traits<Point_3>::Kernel;
  using Segment_3 = typename Kernel::Segment_3;
  std::vector<Segment_3> segments;
  segments.reserve(points.size() - 1);
  for (std::size_t i = 1; i < points.size(); ++i)
    segments.push_back(Segment_3(points[0], points[i]));
  Saver<Kernel> saver;
  saver.export_segments_3(segments, name);
}

template<typename DS, typename CDT>
void dump_cdt(
  const DS& data, const std::size_t sp_idx, const CDT& cdt, std::string file_name) {

  using Point_3       = typename DS::Kernel::Point_3;
  using Vertex_handle = typename CDT::Vertex_handle;

  using Mesh_3 = CGAL::Surface_mesh<Point_3>;
  using VIdx   = typename Mesh_3::Vertex_index;
  using FIdx   = typename Mesh_3::Face_index;
  using UM     = typename Mesh_3::template Property_map<FIdx, unsigned char>;

  Mesh_3 mesh;
  UM red   = mesh.template add_property_map<FIdx, unsigned char>("red"  , 125).first;
  UM green = mesh.template add_property_map<FIdx, unsigned char>("green", 125).first;
  UM blue  = mesh.template add_property_map<FIdx, unsigned char>("blue" , 125).first;

  std::map<Vertex_handle, VIdx> map_v2i;
  for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
    const auto& ipoint = vit->point();
    map_v2i.insert(std::make_pair(
      vit, mesh.add_vertex(data.support_plane(sp_idx).to_3d(ipoint))));
  }

  for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
    std::array<VIdx, 3> vertices;
    for (int i = 0; i < 3; ++i) {
      vertices[i] = map_v2i[fit->vertex(i)];
    }

    const auto face = mesh.add_face(vertices);
    CGAL::Random rand(fit->info().index);
    if (fit->info().index != std::size_t(-1)) {
      red[face]   = (unsigned char)(rand.get_int(32, 192));
      green[face] = (unsigned char)(rand.get_int(32, 192));
      blue[face]  = (unsigned char)(rand.get_int(32, 192));
    }
  }

  file_name += "support-cdt-" + std::to_string(sp_idx) + ".ply";
  std::ofstream out(file_name);
  out.precision(20);
  CGAL::IO::write_PLY(out, mesh);
  out.close();
}

} // namespace internal
} // namespace KSP_3
} // namespace CGAL

#endif // CGAL_KSP_DEBUG_H
