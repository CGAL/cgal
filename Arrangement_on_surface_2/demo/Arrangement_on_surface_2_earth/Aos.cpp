// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "Aos.h"

#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <vector>

#include <qmath.h>
#include <qvector3d.h>

#include <nlohmann/json.hpp>

#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_point_location_result.h>

#include "Aos_defs.h"
#include "arr_print.h"
#include "Tools.h"

using json = nlohmann::ordered_json;

namespace {
  // use this traits every time you construct an arrangement!
  static Geom_traits s_traits;

  // Extended DCEL & Arrangement
  struct Flag {
    bool v;
    Flag() : v{ false } {}
    Flag(bool init) : v{ init } {}
  };

  // EXTENDED AOS for analysing the arrangement
  using Ext_dcel = CGAL::Arr_extended_dcel<Geom_traits, Flag, Flag, Flag>;
  using Ext_topol_traits =
    CGAL::Arr_spherical_topology_traits_2<Geom_traits, Ext_dcel>;
  using Ext_aos = CGAL::Arrangement_on_surface_2<Geom_traits, Ext_topol_traits>;
  using Dir3 = Kernel::Direction_3;
  using Approximate_point_2 = Geom_traits::Approximate_point_2;
  using Approximate_number_type = Geom_traits::Approximate_number_type;
  using Approximate_kernel = Geom_traits::Approximate_kernel;
  using Approximate_Vector_3 = CGAL::Vector_3<Approximate_kernel>;
  using Approximate_Direction_3 = Approximate_kernel::Direction_3;
  using Direction_3 = Kernel::Direction_3;

  //---------------------------------------------------------------------------
  // below are the helper functions used to construct the arcs from KML data
  // TODO: Revisit handling of INNER & OUTER boundaries
  using Curves = std::vector<Curve>;

  // get curves for the given kml placemark
  // NOTE: this is defined here to keep the definitions local to this cpp file
  Curves get_arcs(const Kml::Placemark& placemark) {
    //Geom_traits traits;
    auto ctr_p = s_traits.construct_point_2_object();
    auto ctr_cv = s_traits.construct_curve_2_object();

    std::vector<Curve>  xcvs;
    for (const auto& polygon : placemark.polygons) {
      // colect all rings into a single list (FOR NOW!!!)
      // TO-DO: PROCESS OUTER & INNER BOUNDARIES SEPARATELY!!!
      Kml::LinearRings linear_rings;
      linear_rings.push_back(polygon.outer_boundary);
      for (const auto& inner_boundary : polygon.inner_boundaries)
        linear_rings.push_back(inner_boundary);


      // convert the nodes to points on unit-sphere
      for (const auto& lring : linear_rings) {
        std::vector<Approximate_Vector_3> sphere_points;
        for (const auto& node : lring.nodes) {
          const auto p = node.get_coords_3d();
          Approximate_Vector_3 v(p.x, p.y, p.z);
          sphere_points.push_back(v);
        }

        // add geodesic arcs for the current LinearRing
        std::size_t num_points = sphere_points.size();
        for (std::size_t i = 0; i < num_points - 1; ++i) {
          const auto p1 = sphere_points[i];
          const auto p2 = sphere_points[i + 1];
          xcvs.push_back(ctr_cv(ctr_p(p1.x(), p1.y(), p1.z()),
            ctr_p(p2.x(), p2.y(), p2.z())));
        }
      }
    }

    return xcvs;
  }


  // this one is used by the Aos::check and Aos::ext_check functions
  std::size_t num_counted_nodes = 0;
  std::size_t num_counted_arcs = 0;
  std::size_t num_counted_polygons = 0;
  std::map<Ext_aos::Vertex_handle, Kml::Node>  vertex_node_map;

  template<typename Arr_type>
  Curves  get_arcs(const Kml::Placemarks& placemarks, Arr_type& arr) {
    //Geom_traits traits;
    auto ctr_p = s_traits.construct_point_2_object();
    auto ctr_cv = s_traits.construct_curve_2_object();

    num_counted_nodes = 0;
    num_counted_arcs = 0;
    num_counted_polygons = 0;
    std::vector<Curve>  xcvs;
    for (const auto& pm : placemarks) {
      for (const auto& polygon : pm.polygons) {
        ++num_counted_polygons;

        // colect all rings into a single list (FOR NOW!!!)
        // TO-DO: PROCESS OUTER & INNER BOUNDARIES SEPARATELY!!!
        Kml::LinearRings linear_rings;
        linear_rings.push_back(polygon.outer_boundary);
        for (const auto& inner_boundary : polygon.inner_boundaries)
          linear_rings.push_back(inner_boundary);

        // loop on outer and inner boundaries
        for (const auto& lring : linear_rings) {
          // convert the nodes to points on unit-sphere
          std::vector<Approximate_Vector_3>  sphere_points;
          for (const auto& node : lring.nodes) {
            ++num_counted_nodes;
            const auto p = node.get_coords_3d();
            Approximate_Vector_3  v(p.x, p.y, p.z);
            sphere_points.push_back(v);
            auto vh CGAL_UNUSED = CGAL::insert_point(arr, ctr_p(p.x, p.y, p.z));
            if constexpr (std::is_same<Arr_type, Ext_aos>::value)
            {
              vertex_node_map.insert(std::make_pair(vh, node));
            }
          }

          // add curves
          std::size_t num_points = sphere_points.size();
          for (std::size_t i = 0; i < num_points - 1; ++i) {
            ++num_counted_arcs;
            const auto p1 = sphere_points[i];
            const auto p2 = sphere_points[i + 1];
            auto xcv = ctr_cv(ctr_p(p1.x(), p1.y(), p1.z()),
              ctr_p(p2.x(), p2.y(), p2.z()));
            xcvs.push_back(xcv);
          }
        }
      }
    }
    return xcvs;
  }


  Aos::Approx_arc get_approx_curve(Curve xcv, double error) {
    //Geom_traits traits;
    auto approx = s_traits.approximate_2_object();
    std::vector<QVector3D> approx_curve;
    {
      std::vector<Approximate_point_2> v;
      approx(xcv, error, std::back_insert_iterator(v));

      for (const auto& p : v) {
        const QVector3D arc_point(p.dx(), p.dy(), p.dz());
        approx_curve.push_back(arc_point);
      }
    }

    return approx_curve;
  }
  Aos::Approx_arcs get_approx_curves(std::vector<Curve>& xcvs, double error) {
    Aos::Approx_arcs  approx_curves;
    for (const auto& xcv : xcvs) {
      auto approx_curve = get_approx_curve(xcv, error);
      approx_curves.push_back(std::move(approx_curve));
    }

    return approx_curves;
  }
}

Aos::Approx_arc Aos::get_approx_identification_curve(double error) {
  //Geom_traits traits;
  auto ctr_p = s_traits.construct_point_2_object();
  auto ctr_cv = s_traits.construct_curve_2_object();

  // identification curve (meridian pierced by NEGATIVE Y-AXIS)
  auto xcv = ctr_cv(ctr_p(0, 0, -1), ctr_p(0, 0, 1), Dir3(0, 1, 0));

  auto approx = s_traits.approximate_2_object();
  Approx_arc approx_arc;
  {
    std::vector<Approximate_point_2> v;
    approx(xcv, error, std::back_insert_iterator(v));
    for (const auto& p : v) {
      const QVector3D arc_point(p.dx(), p.dy(), p.dz());
      approx_arc.push_back(arc_point);
    }
  }

  return approx_arc;
}

//! \brief
Aos::Approx_arcs Aos::get_approx_arcs(double error) {
  auto ctr_p = s_traits.construct_point_2_object();
  auto ctr_cv = s_traits.construct_curve_2_object();

  std::vector<Curve>  xcvs;
  xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 1, 0)));
  xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 0, 1)));
  xcvs.push_back(ctr_cv(ctr_p(0, 1, 0), ctr_p(0, 0, 1)));
  //xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 1, 0), Dir3(0, 0, -1)));
  //xcvs.push_back(ctr_cv(Dir3(0, 0, -1)));

  auto approx_arcs = get_approx_curves(xcvs, error);
  return approx_arcs;
}

//! \brief
auto Aos::get_approx_arcs(const Kml::Placemark& placemark, double error) ->
  Approx_arcs {
  auto xcvs = get_arcs(placemark);
  auto arcs = ::get_approx_curves(xcvs, error);
  return arcs;
}

//! \brief
void Aos::check(const Kml::Placemarks& placemarks) {
  //Geom_traits traits;
  Arrangement arr(&s_traits);

  auto xcvs = get_arcs(placemarks, arr);
  std::cout << "-------------------------------\n";
  std::cout << "num arr vertices (before adding arcs) = " <<
    arr.number_of_vertices() << std::endl;
  // add arcs
  for (auto xcv : xcvs)
    CGAL::insert(arr, xcv);

  std::cout << "-------------------------------\n";
  std::cout << "num nodes = " << num_counted_nodes << std::endl;
  std::cout << "num arr vertices = " << arr.number_of_vertices() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num counted arcs = " << num_counted_arcs << std::endl;
  std::cout << "num arr edges = " << arr.number_of_edges() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num polygons = " << num_counted_polygons << std::endl;
  std::cout << "num arr faces = " << arr.number_of_faces() << std::endl;
}

//! \brief
std::vector<QVector3D> Aos::ext_check(const Kml::Placemarks& placemarks) {
  // Construct the arrangement from 12 geodesic arcs.
  Ext_aos arr(&s_traits);

  std::cout << "-------------------------------\n";
  std::cout << "** num arr FACES (before adding arcs) = " <<
    arr.number_of_faces() << std::endl;

  auto xcvs = get_arcs(placemarks, arr);

  // MARK all vertices as true
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    vit->set_data(Flag(true));

  std::cout << "-------------------------------\n";
  std::cout << "num arr vertices (before adding arcs) = " <<
    arr.number_of_vertices() << std::endl;

  // add arcs
  for (auto& xcv : xcvs) CGAL::insert(arr, xcv);

  // extract all vertices that are ADDED when inserting the arcs!
  std::size_t num_created_vertices = 0;
  std::vector<QVector3D> created_vertices;
  auto approx = s_traits.approximate_2_object();
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (vit->data().v == false) {
      std::cout << "-------------------------------------\n";
      std::cout << vit->point() << std::endl;

      if (2 == vit->degree()) {} //continue;

      if (1 == vit->degree()) {
        auto p = vit->point();
        // auto p2 = p.location();
        std::cout << "   deg-1 vertex = " << p << std::endl;
        std::cout << "   deg-1 vertex: " << std::boolalpha
                  << vit->incident_halfedges()->target()->data().v << std::endl;
      }


      ++num_created_vertices;
      auto p = vit->point();
      auto ap = approx(p);
      QVector3D new_vertex(ap.dx(), ap.dy(), ap.dz());
      new_vertex.normalize();
      std::cout << new_vertex << std::endl;
      std::cout << "degree = " << vit->degree() << std::endl;

      created_vertices.push_back(new_vertex);

      // find the arcs that are adjacent to the vertex of degree 4
      if (4 == vit->degree()) {
        std::cout << "**************************\n DEGREE 4 VERTEX: \n";
        const auto first = vit->incident_halfedges();
        auto curr = first;
        do {
          auto tvh = curr->twin()->target();
          //std::cout << std::boolalpha << svh->data().v << " - "
          //          << tvh->data().v << std::endl;
          auto it = vertex_node_map.find(tvh);
          if (it != vertex_node_map.end())
            std::cout << std::setprecision(16) << it->second << std::endl;
          else
            std::cout << "NOT FOUND!!\n";
        } while (++curr != first);
      }

      std::cout << "\n";
    }
  }
  Kml::Node n{ 180.0, -84.71338 };
  std::cout << "Node itself = " << n.get_coords_3d() << std::endl;
  std::cout << "*** num created vertices = " << num_created_vertices
            << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num nodes = " << num_counted_nodes << std::endl;
  std::cout << "num arr vertices = " << arr.number_of_vertices() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num counted arcs = " << num_counted_arcs << std::endl;
  std::cout << "num arr edges = " << arr.number_of_edges() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num polygons = " << num_counted_polygons << std::endl;
  std::cout << "num arr faces = " << arr.number_of_faces() << std::endl;

  return created_vertices;
}

//! \brief
std::vector<QVector3D>  Aos::ext_check_id_based(Kml::Placemarks& placemarks) {
  // Construct the arrangement from 12 geodesic arcs.
  Ext_aos arr(&s_traits);

  auto nodes = Kml::generate_ids(placemarks);
  auto ctr_p = s_traits.construct_point_2_object();
  auto ctr_cv = s_traits.construct_curve_2_object();

  std::size_t num_counted_arcs = 0;
  std::size_t num_counted_polygons = 0;
  //
  std::vector<Point> points;
  std::vector<Ext_aos::Vertex_handle> vertices;
  for (const auto& node : nodes) {
    auto n = node.get_coords_3d();
    auto p = ctr_p(n.x, n.y, n.z);
    auto v = CGAL::insert_point(arr, p);
    points.push_back(p);
    vertices.push_back(v);
    //arr.insert_at_vertices(Segment(p, p), v, v);
  }
  std::cout << "num nodes = " << nodes.size() << std::endl;
  std::cout << "num points = " << points.size() << std::endl;
  // MARK all vertices as true
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    vit->set_data(Flag(true));


  for (auto& placemark : placemarks) {
    for (auto& polygon : placemark.polygons) {
      ++num_counted_polygons;

      // TO DO : ADD the outer boundaries!
      auto& ids = polygon.outer_boundary.ids;
      std::size_t num_nodes = ids.size();
      for (std::size_t i = 0; i < num_nodes - 1; ++i) {
        ++num_counted_arcs;
        const auto nid1 = ids[i];
        const auto nid2 = ids[i + 1];
        auto p1 = points[nid1];
        auto p2 = points[nid2];
        CGAL::insert(arr, ctr_cv(p1, p2));
      }
    }
  }

  std::cout << "-------------------------------\n";
  std::cout << "num arr vertices (before adding arcs) = " <<
    arr.number_of_vertices() << std::endl;

  // extract all vertices that are ADDED when inserting the arcs!
  std::size_t num_created_vertices = 0;
  std::vector<QVector3D> created_vertices;
  auto approx = s_traits.approximate_2_object();
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    // auto& d = vit->data();
    if (vit->data().v == false) {
      std::cout << "-------------------------------------\n";
      std::cout << vit->point() << std::endl;

      if (2 == vit->degree()) {} //continue;

      if (1 == vit->degree()) {
        auto p = vit->point();
        // auto p2 = p.location();
        std::cout << "deg-1 vertex = " << p << std::endl;
        std::cout << "deg-1 vertex: " << std::boolalpha
                  << vit->incident_halfedges()->target()->data().v << std::endl;
      }


      ++num_created_vertices;
      auto p = vit->point();
      auto ap = approx(p);
      QVector3D new_vertex(ap.dx(), ap.dy(), ap.dz());
      new_vertex.normalize();
      std::cout << new_vertex << std::endl;
      std::cout << "degree = " << vit->degree() << std::endl;

      created_vertices.push_back(new_vertex);

      //// find the arcs that are adjacent to this vertex
      //const auto first = vit->incident_halfedges();
      //auto curr = first;
      //do {

      //} while (++curr != first);
      std::cout << std::endl;
    }
  }
  std::cout << "*** num created vertices = " << num_created_vertices
            << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num nodes = " << nodes.size() << std::endl;
  std::cout << "num arr vertices = " << arr.number_of_vertices() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num counted arcs = " << num_counted_arcs << std::endl;
  std::cout << "num arr edges = " << arr.number_of_edges() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num polygons = " << num_counted_polygons << std::endl;
  std::cout << "num arr faces = " << arr.number_of_faces() << std::endl;

  return created_vertices;
}

//! brief
auto Aos::find_new_faces(Kml::Placemarks& placemarks) -> Approx_arcs {
  //Geom_traits traits;
  Ext_aos arr(&s_traits);
  auto ctr_p = s_traits.construct_point_2_object();
  auto ctr_cv = s_traits.construct_curve_2_object();

  auto fh = arr.faces_begin();
  fh->data().v = true;
  std::cout << "num faces = " << arr.number_of_faces() << std::endl;

  auto nodes = Kml::generate_ids(placemarks);

  //-------------------------------------------------------------------------
  // define a set of vertex-handles: use this to check if the face is
  // obtained from the polygon definition, or if it is an additional face
  using Vertex_handle = Ext_aos::Vertex_handle;
  std::map<Vertex_handle, int> vertex_id_map;
  std::set<std::set<int>> all_polygon_node_ids;

  num_counted_nodes = 0;
  num_counted_arcs = 0;
  num_counted_polygons = 0;
  std::vector<Curve>  xcvs;
  for (auto& pm : placemarks) {
    std::cout << pm.name << std::endl;
    for (auto& polygon : pm.polygons) {
      ++num_counted_polygons;

      // colect all rings into a single list (FOR NOW!!!)
      // TO-DO: PROCESS OUTER & INNER BOUNDARIES SEPARATELY!!!
      auto linear_rings = polygon.get_all_boundaries();

      auto* lring = &polygon.outer_boundary;
      {
        std::set<int> polygon_node_ids;

        // convert the nodes to points on unit-sphere
        std::vector<Approximate_Vector_3>  sphere_points;
        //for (const auto& node : lring->nodes)
        //std::cout << "   NUM POLYGON-NODES SIZE = "
        //          << lring->ids.size() << std::endl;
        for (std::size_t i = 0; i < lring->ids.size(); ++i) {
          ++num_counted_nodes;
          const auto id = lring->ids[i];
          const auto& node = lring->nodes[i];
          const auto p = node.get_coords_3d();
          Approximate_Vector_3  v(p.x, p.y, p.z);
          sphere_points.push_back(v);
          auto vh = CGAL::insert_point(arr, ctr_p(p.x, p.y, p.z));
          polygon_node_ids.insert(id);
          vertex_id_map[vh] = id;
          vh->data().v = true;
        }
        //std::cout << "   POLYGON-NODES SET SIZE = "
        //          << polygon_node_ids.size() << std::endl;
        if (lring->ids.size() != (1 + polygon_node_ids.size()))
          std::cout << "*** ASSERTION ERROR!!!!\n";

        all_polygon_node_ids.insert(std::move(polygon_node_ids));

        // add curves
        auto num_points = sphere_points.size();
        for (std::size_t i = 0; i < num_points - 1; ++i) {
          ++num_counted_arcs;
          const auto p1 = sphere_points[i];
          const auto p2 = sphere_points[i + 1];
          auto xcv = ctr_cv(ctr_p(p1.x(), p1.y(), p1.z()),
            ctr_p(p2.x(), p2.y(), p2.z()));
          CGAL::insert(arr, xcv);
        }
      }
    }
  }


  // mark all faces as TRUE (= as existing faces)
  std::size_t num_not_found = 0;
  std::vector<Curve>  new_face_arcs;
  for (auto fh = arr.faces_begin(); fh != arr.faces_end(); ++fh) {
    // skip the spherical face
    std::cout << "num outer_ccbs = " << fh->number_of_outer_ccbs() << std::endl;
    if (fh->number_of_outer_ccbs() == 0) continue;

    // construct the set of all node-ids for the current face
    std::set<int>  face_node_ids_set;
    std::vector<int>  face_node_ids;
    std::vector<Curve>  face_arcs;
    auto first = fh->outer_ccb();
    auto curr = first;
    do {
      auto vh = curr->source();
      // skip if the vertex is due to intersection with the identification curve
      if ((vh->data().v == false) && (vh->degree() == 2)) continue;

      auto id = vertex_id_map[vh];
      face_node_ids_set.insert(id);

      face_arcs.push_back(ctr_cv(curr->source()->point(),
                                 curr->target()->point()));
    } while (++curr != first);
    //std::cout << "counted vertices = " << num_vertices << std::endl;
    //std::cout << "vertices in the set = " << polygon_node_ids.size()
    //          << std::endl;

    auto it = all_polygon_node_ids.find(face_node_ids_set);
    if (it == all_polygon_node_ids.cend()) {
      std::cout << "NOT FOUND!!!\n";
      std::cout << "num nodes = " << face_node_ids_set.size() << std::endl;
      ++num_not_found;
      new_face_arcs.insert(new_face_arcs.end(), face_arcs.begin(),
                           face_arcs.end());
    }
  }
  std::cout << "num not found = " << num_not_found << std::endl;

  auto approx_arcs = get_approx_curves(new_face_arcs, 0.001);
  return approx_arcs;
}

//! \brief
void Aos::save_arr(Kml::Placemarks& placemarks, const std::string& file_name) {
#ifndef USE_EPIC
  //Geom_traits traits;
  Ext_aos arr(&s_traits);
  auto ctr_p = s_traits.construct_point_2_object();
  auto ctr_cv = s_traits.construct_curve_2_object();

  auto fh = arr.faces_begin();
  fh->data().v = true;
  std::cout << "num faces = " << arr.number_of_faces() << std::endl;

  auto nodes = Kml::generate_ids(placemarks);

  //-------------------------------------------------------------------------
  // define a set of vertex-handles: use this to check if the face is
  // obtained from the polygon definition, or if it is an additional face
  using Vertex_handle = Ext_aos::Vertex_handle;
  using Face_handle = Ext_aos::Face_handle;
  std::map<Vertex_handle, int> vertex_id_map;
  std::map<std::set<int>, std::string> all_polygon_node_ids_map;

  // map to associate the created faces with the country names
  // CAUTION: the newly created faces

  num_counted_nodes = 0;
  num_counted_arcs = 0;
  num_counted_polygons = 0;
  std::vector<Curve>  xcvs;
  for (auto& pm : placemarks) {
    std::cout << pm.name << std::endl;
    for (auto& polygon : pm.polygons) {
      ++num_counted_polygons;

      auto* lring = &polygon.outer_boundary;
      {
        // auto num_faces_before = arr.number_of_faces();
        std::set<int> polygon_node_ids;

        // convert the nodes to points on unit-sphere
        std::vector<Approximate_Vector_3> sphere_points;
        //for (const auto& node : lring->nodes)
        //std::cout << "   NUM POLYGON-NODES SIZE = " << lring->ids.size()
        //          << std::endl;
        for (std::size_t i = 0; i < lring->ids.size(); ++i) {
          ++num_counted_nodes;
          const auto id = lring->ids[i];
          const auto& node = lring->nodes[i];
          const auto p = node.get_coords_3d();
          Approximate_Vector_3  v(p.x, p.y, p.z);
          sphere_points.push_back(v);
          auto vh = CGAL::insert_point(arr, ctr_p(p.x, p.y, p.z));
          polygon_node_ids.insert(id);
          vertex_id_map[vh] = id;
          vh->data().v = true;
        }
        //std::cout << "   POLYGON-NODES SET SIZE = "
        //          << polygon_node_ids.size() << std::endl;
        if (lring->ids.size() != (1 + polygon_node_ids.size()))
          std::cout << "*** ASSERTION ERROR!!!!\n";

        all_polygon_node_ids_map.insert(std::make_pair(
          std::move(polygon_node_ids), pm.name));

        // add curves
        std::size_t num_points = sphere_points.size();
        for (std::size_t i = 0; i < num_points - 1; ++i) {
          ++num_counted_arcs;
          const auto p1 = sphere_points[i];
          const auto p2 = sphere_points[i + 1];
          auto xcv = ctr_cv(ctr_p(p1.x(), p1.y(), p1.z()),
            ctr_p(p2.x(), p2.y(), p2.z()));
          CGAL::insert(arr, xcv);
        }
      }
    }
  }

  std::cout << "*** arr.number_of_faces = " << arr.number_of_faces()
            << std::endl;
  std::cout << "*** arr.number_of_halfedges = " << arr.number_of_halfedges()
            << std::endl;
  std::cout << "*** arr.number_of_vertices = " << arr.number_of_vertices()
            << std::endl;

  // DEFINE JSON OBJECT
  json js;
  auto& js_points = js["points"] = json::array();

  ////////////////////////////////////////////////////////////////////////////
  // POINTS
  // define a map from each vertex to its position in the arrangement
  //auto get_num_denum
  using FT = typename Kernel::FT;
  //using json = nlohmann::ordered_json;
  FT ft(0);
  auto ex = ft.exact();
  CGAL::Rational_traits<decltype(ex)> rt;
  typename CGAL::Algebraic_structure_traits<decltype(ex)>::Simplify simplify;

  auto set_num_denum = [&](decltype(ex)& x, json& ratx) {
    simplify(x);
    std::stringstream ss_x_num;
    CGAL::IO::set_ascii_mode(ss_x_num);
    ss_x_num << rt.numerator(x);
    std::string xnum;
    ss_x_num >> xnum;
    ratx["num"] = xnum;

    std::stringstream ss_x_den;
    CGAL::IO::set_ascii_mode(ss_x_den);
    ss_x_den << rt.denominator(x);
    std::string xden;
    ss_x_den >> xden;
    ratx["den"] = xden;
  };

  using Point_ = std::decay_t<decltype(arr.vertices_begin()->point())>;
  std::map<void*, int> point_pos_map;
  std::vector<Point_>  points;
  std::map<Vertex_handle, std::size_t>  vertex_pos_map;
  for (auto vh = arr.vertices_begin(); vh != arr.vertices_end(); ++vh) {
    // add the vertex if not found in the map
    auto it = vertex_pos_map.find(vh);
    if (it == vertex_pos_map.end()) {
      std::size_t new_vh_pos = vertex_pos_map.size();
      vertex_pos_map[vh] = new_vh_pos;

      // write the vertex-data to JSON object
      auto& p = vh->point();
      points.push_back(p);
      point_pos_map[&p] = new_vh_pos;
      auto dx = p.dx().exact();
      auto dy = p.dy().exact();
      auto dz = p.dz().exact();

      json jv;
      jv["location"] = p.location();
      set_num_denum(dx, jv["dx"]);
      set_num_denum(dy, jv["dy"]);
      set_num_denum(dz, jv["dz"]);
      js_points.push_back(std::move(jv));
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // CURVES
  // define a map from each curve to its position in the arrangement
  auto& js_curves = js["curves"] = json::array();
  using Ext_curve = Ext_aos::X_monotone_curve_2;
  std::map<Ext_curve*, std::size_t> curve_pos_map;
  std::size_t num_edges = 0;
  for (auto eh = arr.edges_begin(); eh != arr.edges_end(); ++eh) {
    ++num_edges;
    auto& xcv = eh->curve();
    auto it = curve_pos_map.find(&xcv);
    if (it == curve_pos_map.end()) {
      std::size_t new_xcv_pos = curve_pos_map.size();
      curve_pos_map[&xcv] = new_xcv_pos;

      json je;
      auto& sv = xcv.source();
      auto& tv = xcv.target();
      auto svit = std::find(points.begin(), points.end(), sv);
      auto tvit = std::find(points.begin(), points.end(), tv);
      auto svp = std::distance(points.begin(), svit);
      auto tvp = std::distance(points.begin(), tvit);
      //std::cout << "svp= " << point_pos_map[(void*)&sv] << std::endl;
      je["source"] = svp;
      je["target"] = tvp;
      auto& je_normal = je["normal"];

      // write the vertex-data to JSON object
      auto& n = xcv.normal();
      auto dx = n.dx().exact();
      auto dy = n.dy().exact();
      auto dz = n.dz().exact();
      set_num_denum(dx, je_normal["dx"]);
      set_num_denum(dy, je_normal["dy"]);
      set_num_denum(dz, je_normal["dz"]);

      je["is_vertical"] = xcv.is_vertical();
      je["is_directed_right"] = xcv.is_directed_right();
      je["is_full"] = xcv.is_full();

      js_curves.push_back(std::move(je));
    }
  }
  std::cout << "num edges = " << num_edges << std::endl;
  std::cout << "curve map size = " << curve_pos_map.size() << std::endl;
  std::cout << "js_curves size = " << js_curves.size() << std::endl;

  ////////////////////////////////////////////////////////////////////////////
  // VERTICES
  // there is a one-to-one corresponce between vertices and points
  auto& js_vertices = js["vertices"] = json::array();
  for (auto vh = arr.vertices_begin(); vh != arr.vertices_end(); ++vh) {
    json js_vertex;
    auto vpos = vertex_pos_map[vh];
    js_vertex["point"] = vpos;
    js_vertices.push_back(std::move(js_vertex));
  }

  ////////////////////////////////////////////////////////////////////////////
  // EDGES
  num_edges = 0;
  auto& js_edges = js["edges"] = json::array();
  std::map<void*, std::size_t> halfedge_pos_map;
  for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    auto& edge = *eit;
    auto& xcv = edge.curve();
    auto xcvp = curve_pos_map[&xcv];
    json js_edge;
    js_edge["curve"] = xcvp;
    js_edge["direction"] = edge.direction();
    js_edge["source"] = vertex_pos_map[edge.source()];
    js_edge["target"] = vertex_pos_map[edge.target()];
    js_edges.push_back(std::move(js_edge));
    ++num_edges;

    // add the halfedge indices to the map
    std::size_t new_halfedge_index = halfedge_pos_map.size();
    auto& twin = *edge.twin();
    halfedge_pos_map[&edge] = new_halfedge_index;
    halfedge_pos_map[&twin] = new_halfedge_index + 1;
  }
  std::cout << "EDGE CHECKS:\n";
  std::cout << "  *** num edges = " << num_edges << std::endl;
  std::cout << " *** js_edges size = " << js_edges.size() << std::endl;


  ////////////////////////////////////////////////////////////////////////////
  // FACES
  // CONDITION DATA: Caspian Sea needs to be defined
  num_edges = 0;
  std::size_t num_not_found = 0;
  std::map<Face_handle, std::string>  face_name_map;
  for (auto fh = arr.faces_begin(); fh != arr.faces_end(); ++fh) {
    // skip the spherical face
    std::cout << "num outer_ccbs = " << fh->number_of_outer_ccbs() << std::endl;
    if (fh->number_of_outer_ccbs() == 0) continue;

    // construct the set of all node-ids for the current face
    std::set<int> face_node_ids_set;
    std::vector<int> face_node_ids;
    auto first = fh->outer_ccb();
    auto curr = first;
    do {
      auto vh = curr->source();
      // skip if the vertex is due to intersection with the identification curve
      if ((vh->data().v == false) && (vh->degree() == 2)) continue;

      auto id = vertex_id_map[vh];
      face_node_ids_set.insert(id);
    } while (++curr != first);

    std::string name;
    auto it = all_polygon_node_ids_map.find(face_node_ids_set);
    if (it == all_polygon_node_ids_map.cend())
    {
      std::cout << "NOT FOUND!!!\n";
      std::cout << "num nodes = " << face_node_ids_set.size() << std::endl;
      ++num_not_found;
      name = "Caspian Sea";
    }
    else {
      // ++num_found;
      name = it->second;
    }
    face_name_map[fh] = name;
  }
  std::cout << "num not found = " << num_not_found << std::endl;


  // RECORD FACES
  json& js_faces = js["faces"] = json::array();
  auto get_ccb_json = [&](Ext_aos::Ccb_halfedge_circulator first) {
    json js_edges;
    auto& ccb_edge_indices = js_edges["halfedges"] = json::array();
    auto curr = first;
    do {
      auto& he = *curr;
      auto it = halfedge_pos_map.find(&he);
      if (it == halfedge_pos_map.end())
      {
        std::cout << "ASSERTION ERROR!!!" << std::endl;
      }

      auto edge_pos = it->second;
      ccb_edge_indices.push_back(edge_pos);
    } while (++curr != first);

    return js_edges;
  };

  std::size_t total_num_half_edges = 0;
  for (auto fh = arr.faces_begin(); fh != arr.faces_end(); ++fh) {
    //// skip the spherical face
    //if (fh->number_of_outer_ccbs() == 0)
    //  continue;

    // json object for the current face
    json js_face;
    auto face_name = face_name_map[fh];
    js_face["name"] = face_name;
    js_face["is_unbounded"] = false;
    js_face["is_valid"] = true;

    // at this point we are sure that we have at least 1 outer-ccb
    auto& js_outer_ccbs = js_face["outer_ccbs"] = json::array();
    for (auto ccb = fh->outer_ccbs_begin(); ccb != fh->outer_ccbs_end(); ++ccb)
    {
      auto js_ccb = get_ccb_json(*ccb);
      js_outer_ccbs.push_back(std::move(js_ccb));
    }

    // INNER CCBS
    if (fh->number_of_inner_ccbs() > 0) {
      auto& js_inner_ccbs = js_face["inner_ccbs"] = json::array();
      for (auto ccb = fh->inner_ccbs_begin(); ccb != fh->inner_ccbs_end();
           ++ccb) {
        auto js_ccb = get_ccb_json(*ccb);
        js_inner_ccbs.push_back(std::move(js_ccb));
      }
    }

    js_faces.push_back(std::move(js_face));
  }
  std::cout << "total num half-edges = " << total_num_half_edges << std::endl;

  // save the arrangement
  std::ofstream ofile(file_name);
  ofile << js.dump(2);
  ofile.close();
#endif
}

namespace {
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
    Illegal_input(Error_id /* err */, const std::string& msg,
      const std::string& filename) :
      std::logic_error(std::string(msg).append(" (").append(filename).
        append(")!"))
    {}

    Illegal_input(Error_id /* err */, const std::string& msg) :
      std::logic_error(std::string(msg).append("!"))
    {}
  };

  struct Input_file_missing_error : public std::logic_error {
    Input_file_missing_error(std::string& str) : std::logic_error(str) {}
  };

  //
  struct country {
    //! Constructor
    country(std::string& name) : m_name(std::move(name)) {}

    std::string m_name;
  };

  /*! Read a json file.
   */
  bool read_json(const std::string& filename, nlohmann::json& data) {
    using json = nlohmann::json;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
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

  /*!
   */
  template <typename FT>
  FT to_ft(const nlohmann::json& js_ft) {
    using Exact_type = typename FT::Exact_type;
    const std::string& js_num = js_ft["num"];
    const std::string& js_den = js_ft["den"];
    std::string str = js_num + "/" + js_den;
    Exact_type eft(str);
    return FT(eft);
  }

  template <typename Kernel_, typename Arrangement_>
  bool read_arrangement(const std::string& filename, Arrangement_& arr) {
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

    std::cout << "# points: " << num_points << std::endl;
    std::cout << "# curves: " << num_curves << std::endl;
    std::cout << "# vertices: " << num_vertices << std::endl;
    std::cout << "# halfedges: " << num_halfedges << std::endl;
    std::cout << "# faces: " << num_faces << std::endl;

    using Point = typename Arrangement::Point_2;
    using X_monotone_curve = typename Arrangement::X_monotone_curve_2;
    using FT = typename Kernel::FT;
    // using Exact_type = typename FT::Exact_type;

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
      auto direction = js_edge["direction"];
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
    const bool is_unbounded(false);
    const bool is_valid(true);
    for (const auto& js_face : js_faces) {
      DFace* new_f = arr_access.new_face();

      new_f->set_unbounded(is_unbounded);
      new_f->set_fictitious(!is_valid);
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
            prev_he->set_next(curr_he); // connect
            curr_he->set_outer_ccb(new_occb);// set the CCB
            prev_he = curr_he;
          }
          prev_he->set_next(first_he);// close the loop
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
            prev_he->set_next(curr_he);        // connect
            curr_he->set_inner_ccb(new_iccb);  // set the CCB
            prev_he = curr_he;
          }
          prev_he->set_next(first_he);         // close the loop
          new_f->add_inner_ccb(new_iccb, first_he);
        }
      }
    }

    return true;
  }
}

namespace {

  /*!
   */
  template<typename T>
  Aos::Arr_handle get_handle(T* arr) {
    return Aos::Arr_handle(arr, [](void* ptr) {
      // std::cout << "DELETING THE ARRANGEMENT WITH SHARED_PTR DELETER!!\n";
      delete reinterpret_cast<T*>(ptr);
    });
  }
}

/*!
 */
Aos::Arr_handle Aos::load_arr(const std::string& file_name) {
  auto* arr = new Countries_arr(&s_traits);
  auto arrh = get_handle(arr);
  auto rc = read_arrangement<Kernel>(file_name, *arr);
  if (! rc) return nullptr;
  return arrh;
}

Aos::Arr_handle  Aos::construct(Kml::Placemarks& placemarks) {
  auto* arr = new Arrangement(&s_traits);
  auto arrh = get_handle(arr);
  auto xcvs = get_arcs(placemarks, *arr);
  for (auto& xcv : xcvs) CGAL::insert(*arr, xcv);
  return arrh;
}

Aos::Country_color_map  Aos::get_color_mapping(Arr_handle arrh) {
  auto& arr = *reinterpret_cast<Countries_arr*>(arrh.get());

  // group the faces by their country name,
  std::vector<std::string> all_countries;
  using Face_ = Countries_arr::Face_handle::value_type;
  std::map<std::string, std::vector<Face_*>>  country_faces_map;
  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    auto& face = *fit;
    const auto& country_name = fit->data();
    // skipping spherical-face
    if (country_name.empty()) continue;
    country_faces_map[country_name].push_back(&face);
    all_countries.push_back(country_name);
  }

  // prepare a map of neighboring countries
  std::map<std::string, std::set<std::string>> country_neighbors_map;
  for (auto& [country_name, faces] : country_faces_map) {
    // loop on all of the faces of the current country
    for (auto* face : faces) {
      auto first = face->outer_ccb();
      auto curr = first;
      do {
        const auto& neighbor_country_name = curr->twin()->face()->data();

        // skip the spherical face
        if (neighbor_country_name.empty()) continue;

        country_neighbors_map[country_name].insert(neighbor_country_name);
      } while (++curr != first);
    }
  }

  // find a color index for each country by looking at its neighbors
  Country_color_map result;
  for (const auto& country_name : all_countries) {
    // first: find a free color index
    bool color_used[5] = { false, false, false, false, false };
    auto& neighbor_set = country_neighbors_map[country_name];
    for (auto& neighbor : neighbor_set) {
      auto it = result.find(neighbor);
      // if there is a country in the map, then it must have been assigned one!
      if (it != result.end()) {
        auto used_color_index = it->second;
        color_used[used_color_index] = true;
      }
    }

    // find the first color index not used
    bool found = false;
    for (std::size_t i = 0; i < 5; ++i) {
      if (! color_used[i]) {
        found = true;
        result[country_name] = i;
      }
    }
    // assertion check!!!
    if(! found) std::cout << "*** ASSERTION ERROR: NO INDEX FOUND!!!\n";
  }

  return result;
}

std::string Aos::locate_country(Arr_handle arrh, const QVector3D& point) {
  using Naive_pl = CGAL::Arr_naive_point_location<Countries_arr>;

  auto& arr = *reinterpret_cast<Countries_arr*>(arrh.get());
  auto ctr_point = arr.geometry_traits()->construct_point_2_object();
  auto query_point = ctr_point(point.x(), point.y(), point.z());

  Naive_pl npl(arr);
  auto obj = npl.locate(query_point);

  using Arrangement_2 = Countries_arr;
  using Vertex_const_handle = typename Arrangement_2::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement_2::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement_2::Face_const_handle;
  std::string country_name = "";
  if (auto f = std::get_if<Face_const_handle>(&obj)) { // located inside a face
    //std::cout << "*** QUERY: FACE\n";
    country_name = f->ptr()->data();
  }
  else if (auto e = std::get_if<Halfedge_const_handle>(&obj)) {
    // located on an edge: return one of the incident face arbitrarily
    //std::cout << "*** QUERY: EDGE\n";
    country_name = (*e)->face()->data();
  }
  else if (auto v = std::get_if<Vertex_const_handle>(&obj)) {
    // located on a vertex
    if ((*v)->is_isolated()) {
      //std::cout << "*** QUERY: ISOLATED VERTEX\n";
      country_name = (*v)->face()->data();
    }
    else {
      //std::cout << "*** QUERY: VERTEX\n";
      country_name = (*v)->incident_halfedges()->face()->data();
    }
  }
  else CGAL_error_msg("Invalid object.");

  return country_name;
}

Aos::Approx_arcs
Aos::get_approx_arcs_from_faces_edges(Arr_handle arrh, double error) {
  auto& arr = *reinterpret_cast<Countries_arr*>(arrh.get());
  auto ctr_cv = s_traits.construct_curve_2_object();
  Curves xcvs;
  for (auto eit = arr.halfedges_begin(); eit != arr.halfedges_end(); ++eit) {
    auto& s = eit->curve();
    xcvs.push_back( ctr_cv(s.source(), s.target()) );
  }

  auto arcs = ::get_approx_curves(xcvs, error);
  return arcs;
}
