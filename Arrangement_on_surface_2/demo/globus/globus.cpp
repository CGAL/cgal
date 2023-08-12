// Copyright (c) 2022, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel <efifogel@gmail.com>

#include <string>
#include <vector>

#include <QApplication>

// #include <FileGDBAPI.h>
#include <nlohmann/json.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/draw_arrangement_2.h>
#include <CGAL/Arr_accessor.h>

#include "Globus_window.h"

enum Error_id {
  FILE_NOT_FOUND,
  ILLEGAL_EXTENSION,
  UNABLE_TO_OPEN,
  FILE_IS_EMPTY,
  INVALID_INITIAL_POLYGON,
  UNSUPPORTED,
  INVALID_OUTPUT_FILE,
  ILLEGAL_SOLUTION_FILE
};

struct Illegal_input : public std::logic_error {
  Illegal_input(Error_id /* err */, const std::string &msg,
    const std::string &filename) :
    std::logic_error(std::string(msg).append(" (").append(filename).
      append(")!"))
  {}

  Illegal_input(Error_id /* err */, const std::string &msg) :
    std::logic_error(std::string(msg).append("!"))
  {}
};

struct Input_file_missing_error : public std::logic_error {
  Input_file_missing_error(std::string &str) : std::logic_error(str) {}
};

//
struct country {
  //! Constructor
  country(std::string& name) :
    m_name(std::move(name))
  {}

  std::string m_name;
};

// namespace FGA = FileGDBAPI;

#if 0
/*!
 */
int load_countries(std::vector<country>& countries) {
  // Open the geodatabase and the table.
  fgdbError hr;
  FGA::Geodatabase geodatabase;
  hr = FGA::OpenGeodatabase(L"/home/efif/tmp/esri/19113835-45b5-40b9-989b-fd04aad324da.gdb", geodatabase);
  if (hr != S_OK) return hr;

  // std::vector<std::wstring> types;
  // hr = geodatabase.GetDatasetTypes(types);
  // if (hr != S_OK) return hr;
  // std::cout << "Types:\n";
  // for (const auto& type : types)
  //   std::wcout << type << std::endl;

  FGA::Table table;
  hr = geodatabase.OpenTable(L"World_Countries", table);
  if (hr != S_OK) {
    std::cerr << "Cannot open table!\n";
    return hr;
  }

  // std::string doc;
  // hr = table.GetDocumentation(doc);
  // if (hr != S_OK) return hr;
  // std::cout << "Table Information:\n";
  // std::cout << doc << std::endl;

  // FGA::fieldDefs field_defs;
  // hr = table.GetFields(fields_defs);
  // if (hr != S_OK) return hr;
  // std::cout << "Fields:\n";
  // for (const auto& field : Fields)
  //   std::wcout << field << std::endl;
  FGA::FieldInfo field_info;
  hr = table.GetFieldInformation(field_info);
  if (hr != S_OK) return hr;
  int count;
  hr = field_info.GetFieldCount(count);
  std::cout << "Number of fileds: " << count << std::endl;
  for (auto i = 0; i < count; ++i) {
    std::wstring name;
    hr = field_info.GetFieldName(i, name);
    if (hr != S_OK) return hr;
    std::wcout << "Name: " << name << std::endl;
  }

  FGA::Envelope envelope;
  FGA::EnumRows rows;
  hr = table.Search(L"SHAPE, COUNTRY", L"", envelope, true, rows);
  if (hr != S_OK) {
    std::cerr << "Cannot find rows!\n";
    return hr;
  }

  countries.clear();
  FGA::Row row;
  while (rows.Next(row) == S_OK) {
    std::wstring name;
    row.GetString(L"COUNTRY", name);
    std::string simple_name;
    simple_name.assign(name.begin(), name.end());
    countries.emplace_back(simple_name);
  }

  // Close the table
  hr = geodatabase.CloseTable(table);
  if (hr != S_OK) {
    std::wcout << "An error occurred while closing Cities." << endl;
    std::wcout << "Error code: " << hr << endl;
    return -1;
  }

  // Close the geodatabase
  hr = FGA::CloseGeodatabase(geodatabase);
  if (hr != S_OK) {
    std::wcout << "An error occurred while closing the geodatabase." << endl;
    std::wcout << "Error code: " << hr << endl;
    return -1;
  }

  return 0;
}
#endif

/*! Read a json file.
 */
bool read_json(const std::string& filename, nlohmann::json& data) {
  using json = nlohmann::json;
  std::ifstream infile(filename);
  if (! infile.is_open()) {
    throw Illegal_input(UNABLE_TO_OPEN, "Cannot open file", filename);
    return false;
  }
  data = json::parse(infile);
  infile.close();
  if (data.empty()) {
    throw Illegal_input(FILE_IS_EMPTY, "File is empty", filename);
    return false;
  }
  return true;
}

template <typename FT>
FT to_ft(const nlohmann::json& js_ft) {
  using Exact_type = typename FT::Exact_type;
  const std::string& js_num = js_ft["num"];
  const std::string& js_den = js_ft["den"];
  std::string str = js_num + "/" + js_den;
  Exact_type eft(str);
  return FT(eft);
}

template <typename Arrangement_, typename Kernel_>
bool read_arrangement(const std::string& filename, Arrangement_& arr,
                      const Kernel_& kernel) {
  using Arrangement = Arrangement_;
  using Kernel = Kernel_;

  using json = nlohmann::json;
  json data;
  auto rc = read_json(filename, data);
  if (! rc) return false;

  // points
  auto it = data.find("points");
  if (it == data.end()) {
    std::cerr << "The points item is missing " << " (" << filename << ")\n";
    return false;
  }
  const auto& js_points = it.value();

  // curves
  it = data.find("curves");
  if (it == data.end()) {
    std::cerr << "The curves item is missing " << " (" << filename << ")\n";
    return false;
  }
  const auto& js_curves = it.value();

  // vertices
  it = data.find("vertices");
  if (it == data.end()) {
    std::cerr << "The vertices item is missing " << " (" << filename << ")\n";
    return false;
  }
  const auto& js_vertices = it.value();

  // edges
  it = data.find("edges");
  if (it == data.end()) {
    std::cerr << "The edges item is missing " << " (" << filename << ")\n";
    return false;
  }
  const auto& js_edges = it.value();

  // faces
  it = data.find("faces");
  if (it == data.end()) {
    std::cerr << "The faces item is missing " << " (" << filename << ")\n";
    return false;
  }
  const auto& js_faces = it.value();

  const std::size_t num_points = js_points.size();
  const std::size_t num_curves = js_curves.size();
  const std::size_t num_vertices = js_vertices.size();
  const std::size_t num_edges = js_edges.size();
  const std::size_t num_faces = js_faces.size();
  const std::size_t num_halfedges = num_edges * 2;

  if (num_points < num_vertices) {
    std::cerr << "The no. of points (" << num_points
              << ") is smaller than the no. of vertices (" << num_vertices
              << ")\n";
    return false;
  }

  if (num_curves < num_edges) {
    std::cerr << "The no. of curves (" << num_curves
              << ") is smaller than the no. of edge (" << num_edges << ")\n";
    return false;
  }

  std::cout << "# points: " << num_points<< std::endl;
  std::cout << "# curves: " <<  num_curves<< std::endl;
  std::cout << "# vertices: " << num_vertices << std::endl;
  std::cout << "# halfedges: " << num_halfedges << std::endl;
  std::cout << "# faces: " << num_faces << std::endl;

  using Point = typename Arrangement::Point_2;
  using X_monotone_curve = typename Arrangement::X_monotone_curve_2;
  using FT = typename Kernel::FT;
  using Exact_type = typename FT::Exact_type;

  std::vector<Point> points;
  points.reserve(num_points);
  for (const auto& js_pnt : js_points) {
    using Direction_3 = typename Kernel::Direction_3;
    using Location = typename Point::Location_type;
    auto location = static_cast<Location>(js_pnt["location"]);
    auto dx = to_ft<FT>(js_pnt["dx"]);
    auto dy = to_ft<FT>(js_pnt["dy"]);
    auto dz = to_ft<FT>(js_pnt["dz"]);
    Direction_3 dir(dx, dy, dz);
    Point pnt(dir, location);
    points.push_back(pnt);
  }

  std::vector<X_monotone_curve> xcurves;
  xcurves.reserve(num_curves);
  for (const auto& js_xcv : js_curves) {
    using Direction_3 = typename Kernel::Direction_3;
    std::size_t src_id = js_xcv["source"];
    std::size_t trg_id = js_xcv["target"];
    const auto& js_normal = js_xcv["normal"];
    auto dx = to_ft<FT>(js_normal["dx"]);
    auto dy = to_ft<FT>(js_normal["dy"]);
    auto dz = to_ft<FT>(js_normal["dz"]);
    Direction_3 normal(dx, dy, dz);
    bool is_vert = js_xcv["is_vertical"];
    bool is_directed_right = js_xcv["is_directed_right"];
    bool is_full = js_xcv["is_full"];
    const auto& src = points[src_id];
    const auto& trg = points[trg_id];
    X_monotone_curve xcv(src, trg, normal, is_vert, is_directed_right, is_full);
    xcurves.push_back(xcv);
  }

  using Arr_accessor = CGAL::Arr_accessor<Arrangement>;
  Arr_accessor arr_access(arr);
  arr_access.clear_all();

  // Vertices
  using DVertex = typename Arr_accessor::Dcel_vertex;
  std::vector<DVertex*> vertices(num_vertices);
  size_t k = 0;
  for (const auto& js_vertex : js_vertices) {
    std::size_t point_id = js_vertex["point"];
    const auto& point = points[point_id];
    CGAL::Arr_parameter_space ps_x, ps_y;
    switch (point.location()) {
     case Point::NO_BOUNDARY_LOC: ps_x = ps_y = CGAL::INTERIOR; break;
     case Point::MIN_BOUNDARY_LOC:
      ps_x = CGAL::INTERIOR;
      ps_y = CGAL::ARR_BOTTOM_BOUNDARY;
      break;
     case Point::MID_BOUNDARY_LOC:
      ps_x = CGAL::LEFT_BOUNDARY;
      ps_y = CGAL::INTERIOR;
      break;
     case Point::MAX_BOUNDARY_LOC:
      ps_x = CGAL::INTERIOR;
      ps_y = CGAL::ARR_TOP_BOUNDARY;
      break;
    }
    vertices[k++] = arr_access.new_vertex(&point, ps_x, ps_y);
  }

  // Halfedges
  using DHalfedge = typename Arr_accessor::Dcel_halfedge;
  std::vector<DHalfedge*> halfedges(num_halfedges);
  k = 0;
  for (const auto& js_edge : js_edges) {
    std::size_t source_id = js_edge["source"];
    std::size_t target_id = js_edge["target"];
    std::size_t curve_id = js_edge["curve"];
    int direction = js_edge["direction"];
    DVertex* src_v = vertices[source_id];
    DVertex* trg_v = vertices[target_id];
    const auto& curve = xcurves[curve_id];
    DHalfedge* new_he = arr_access.new_edge(&curve);
    trg_v->set_halfedge(new_he);
    new_he->set_vertex(trg_v);
    src_v->set_halfedge(new_he->opposite());
    new_he->opposite()->set_vertex(src_v);
    new_he->set_direction(static_cast<CGAL::Arr_halfedge_direction>(direction));
    halfedges[k++] = new_he;
    halfedges[k++] = new_he->opposite();
  }

  // Faces
  using DFace = typename Arr_accessor::Dcel_face;
  using DOuter_ccb = typename Arr_accessor::Dcel_outer_ccb;
  using DInner_ccb = typename Arr_accessor::Dcel_inner_ccb;
  using DIso_vert = typename Arr_accessor::Dcel_isolated_vertex;
  const bool is_unbounded(false);
  const bool is_valid(true);
  for (const auto& js_face : js_faces) {
    DFace* new_f = arr_access.new_face();

    new_f->set_unbounded(is_unbounded);
    new_f->set_fictitious(! is_valid);
    new_f->set_data(js_face["name"]);
    // Read the outer CCBs of the face.
    auto oit = js_face.find("outer_ccbs");
    if (oit != js_face.end()) {
      const auto& js_outer_ccbs = *oit;
      for (const auto& js_ccb : js_outer_ccbs) {
        // Allocate a new outer CCB record and set its incident face.
        auto* new_occb = arr_access.new_outer_ccb();
        new_occb->set_face(new_f);

        // Read the current outer CCB.
        auto bit = js_ccb.find("halfedges");
        if (bit == js_ccb.end()) {
          std::cerr << "The halfedges item is missing " << " (" << filename
                    << ")\n";
          return false;
        }

        const auto& js_halfedges = *bit;
        auto hit = js_halfedges.begin();
        std::size_t first_idx = *hit;
        DHalfedge* first_he = halfedges[first_idx];
        first_he->set_outer_ccb(new_occb);
        DHalfedge* prev_he = first_he;
        for (++hit; hit != js_halfedges.end(); ++hit) {
          std::size_t curr_idx = *hit;
          auto curr_he = halfedges[curr_idx];
          prev_he->set_next(curr_he);		// connect
          curr_he->set_outer_ccb(new_occb);	// set the CCB
          prev_he = curr_he;
        }
        prev_he->set_next(first_he);		// close the loop
        new_f->add_outer_ccb(new_occb, first_he);
      }
    }

    // Read the inner CCBs of the face.
    auto iit = js_face.find("inner_ccbs");
    if (iit != js_face.end()) {
      const auto& js_inner_ccbs = *iit;
      for (const auto& js_ccb : js_inner_ccbs) {
        // Allocate a new inner CCB record and set its incident face.
        auto* new_iccb = arr_access.new_inner_ccb();
        new_iccb->set_face(new_f);

        // Read the current inner CCB.
        auto bit = js_ccb.find("halfedges");
        if (bit == js_ccb.end()) {
          std::cerr << "The halfedges item is missing " << " (" << filename
                    << ")\n";
          return false;
        }

        const auto& js_halfedges = *bit;
        auto hit = js_halfedges.begin();
        std::size_t first_idx = *hit;
        DHalfedge* first_he = halfedges[first_idx];
        first_he->set_inner_ccb(new_iccb);
        DHalfedge* prev_he = first_he;
        for (++hit; hit != js_halfedges.end(); ++hit) {
          std::size_t curr_idx = *hit;
          auto curr_he = halfedges[curr_idx];
          prev_he->set_next(curr_he);		// connect
          curr_he->set_inner_ccb(new_iccb);	// set the CCB
          prev_he = curr_he;
        }
        prev_he->set_next(first_he);		// close the loop
        new_f->add_inner_ccb(new_iccb, first_he);
      }
    }

    // // Read the isolated vertices inside the face.
    // Size n_isolated_vertices =
    //   formatter.read_size("number_of_isolated_vertices");
    // if (n_isolated_vertices) {
    //   formatter.read_isolated_vertices_begin();
    //   Size k;
    //   for (k = 0; k < n_isolated_vertices; k++) {
    //     // Allocate a new isolated vertex record and set its incident face.
    //     DIso_vert* new_iso_vert = arr_access.new_isolated_vertex();
    //     new_iso_vert->set_face(new_f);
    //     // Read the current isolated vertex.
    //     std::size_t v_idx = formatter.read_vertex_index();
    //     DVertex* iso_v = m_vertices[v_idx];
    //     iso_v->set_isolated_vertex(new_iso_vert);
    //     new_f->add_isolated_vertex(new_iso_vert, iso_v);
    //   }
    // }
  }

  return true;
}

//
int main(int argc, char* argv[]) {
  const char* filename = (argc > 1) ? argv[1] : "arr.json";

  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Geom_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
  using Point = Geom_traits::Point_2;
  using X_monotone_curve = Geom_traits::X_monotone_curve_2;
  using Dcel = CGAL::Arr_face_extended_dcel<Geom_traits, std::string>;
  using Topol_traits = CGAL::Arr_spherical_topology_traits_2<Geom_traits, Dcel>;
  using Arrangement = CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>;

  Kernel kernel;
  Geom_traits traits;
  Arrangement arr(&traits);;
  auto rc = read_arrangement(filename, arr, kernel);
  if (! rc) {
    std::cerr << "Failed to load database!\n";
    return -1;
  }
  std::cout << arr << std::endl;
  CGAL::draw(arr);

  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("CGAL");
  QCoreApplication::setApplicationName("globus");

  // Import resources from libCGAL (Qt5).
  CGAL_QT_INIT_RESOURCES;

  Globus_window window;
  window.show();
  return app.exec();
}
