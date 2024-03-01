// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <vector>
#include <unordered_map>

#include <boost/property_map/property_map.hpp>

#include <qmath.h>
#include <qvector3d.h>

#include <nlohmann/json.hpp>

// Includes for Constrained Delaunay Triangulation
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Polygon_2.h>
//#include <CGAL/Projection_traits_3.h>

#include "Aos_defs.h"
#include "Aos_triangulator.h"
//#include "arr_print.h"
//#include "Tools.h"

using json = nlohmann::ordered_json;

namespace {

  // use this traits everytime you construct an arrangment!
  static Geom_traits s_traits;

  using Dir3 = Kernel::Direction_3;
  std::ostream& operator << (std::ostream& os, const Dir3& d) {
    os << d.dx() << ", " << d.dy() << ", " << d.dz();
    return os;
  }

  using Approximate_point_2 = Geom_traits::Approximate_point_2;
  std::ostream& operator << (std::ostream& os, const Approximate_point_2& d) {
    os << d.dx() << ", " << d.dy() << ", " << d.dz();
    return os;
  }

  using Approximate_number_type = Geom_traits::Approximate_number_type;
  using Approximate_kernel = Geom_traits::Approximate_kernel;
  using Approximate_Vector_3 = CGAL::Vector_3<Approximate_kernel>;
  using Approximate_Direction_3 = Approximate_kernel::Direction_3;
  using Direction_3 = Kernel::Direction_3;
}

//! \brief
std::vector<QVector3D> Aos_triangulator::get_all(Aos::Arr_handle arrh) {
  using K     = CGAL::Exact_predicates_inexact_constructions_kernel;
  //using K     = CGAL::Projection_traits_3<K_epic>;
  using Vb    = CGAL::Triangulation_vertex_base_2<K>;
  using Fb    = CGAL::Constrained_triangulation_face_base_2<K>;
  using TDS   = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Itag  = CGAL::Exact_predicates_tag;
  using CDT   = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>;
  using Face_handle = CDT::Face_handle;
  using Point       = CDT::Point;
  using Polygon_2   = CGAL::Polygon_2<K>;

  auto& arr = *reinterpret_cast<Arrangement*>(arrh.get());

  //Geom_traits traits;
  auto approx = s_traits.approximate_2_object();


  std::vector<std::vector<QVector3D>> all_faces;
  // loop on all faces of the arrangement
  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    // skip any face with no OUTER-CCB
    if (0 == fit->number_of_outer_ccbs())
      continue;

    // COMPUTE THE CENTROID OF ALL FACE-POINTS
    std::vector<QVector3D> face_points;

    // loop on the egdes of the current outer-ccb
    auto first = fit->outer_ccb();
    auto curr = first;
    do {
      auto ap = approx(curr->source()->point());
      QVector3D p(ap.dx(), ap.dy(), ap.dz());
      p.normalize();
      face_points.push_back(p);
    } while (++curr != first);

    all_faces.push_back(std::move(face_points));
  }

  // RESULTING TRIANGLE POINTS (every 3 point => triangle)
  std::vector<QVector3D>  triangles;

  std::cout << "triangulating individual faces\n";

  // loop on all approximated faces
  for (auto& face_points : all_faces) {
    //std::cout << "num face points = " << face_points.size() << std::endl;
    // no need to triangulate if the number of points is 3
    if (face_points.size() == 3) {
      triangles.insert(triangles.end(), face_points.begin(), face_points.end());
      continue;
    }

    // find the centroid of all face-points
    QVector3D centroid(0, 0, 0);
    for (const auto& fp : face_points) centroid += fp;
    centroid /= face_points.size();
    centroid.normalize();
    auto normal = centroid;

    K::Point_3  plane_origin(centroid.x(), centroid.y(), centroid.z());
    K::Vector_3 plane_normal(normal.x(), normal.y(), normal.z());
    K::Plane_3 plane(plane_origin, plane_normal);

    Polygon_2 polygon;

    // project all points onto the plane
    K::Point_3 origin(0, 0, 0);
    for (const auto& fp : face_points) {
      // define a ray through the origin and the current point
      K::Point_3 current_point(fp.x(), fp.y(), fp.z());
      K::Ray_3 ray(origin, current_point);

      auto intersection = CGAL::intersection(plane, ray);
      if (!intersection.has_value())
        std::cout << "INTERSECTION ASSERTION ERROR!!!\n";
      auto ip = std::get<K::Point_3>(intersection.value());
      auto ip2 = plane.to_2d(ip);

      // add this to the polygon constraint
      polygon.push_back(ip2);
    }

    CDT cdt;
    cdt.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(),
                          true);

    std::unordered_map<Face_handle, bool> in_domain_map;
    boost::associative_property_map< std::unordered_map<Face_handle, bool>>
      in_domain(in_domain_map);

    //Mark facets that are inside the domain bounded by the polygon
    CGAL::mark_domain_in_triangulation(cdt, in_domain);

    // loop on all the triangles ("faces" in triangulation doc)
    for (Face_handle f : cdt.finite_face_handles()) {
      // if the current triangles is not inside the polygon -> skip it
      if (false == get(in_domain, f)) continue;

      for (int i = 0; i < 3; ++i) {
        auto tp = f->vertex(i)->point();
        auto tp3 = plane.to_3d(tp);
        QVector3D p3(tp3.x(), tp3.y(), tp3.z());
        p3.normalize();
        triangles.push_back(p3);
      }
    }
  }

  return triangles;
}

//! \brief
Country_triangles_map Aos_triangulator::get_by_country(Aos::Arr_handle arrh) {
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  //using K     = CGAL::Projection_traits_3<K_epic>;
  using Vb = CGAL::Triangulation_vertex_base_2<K>;
  using Fb = CGAL::Constrained_triangulation_face_base_2<K>;
  using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Itag = CGAL::Exact_predicates_tag;
  using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>;
  using Face_handle = CDT::Face_handle;
  using Point = CDT::Point;
  using Polygon_2 = CGAL::Polygon_2<K>;

  auto& arr = *reinterpret_cast<Countries_arr*>(arrh.get());
  auto approx = s_traits.approximate_2_object();

  // group the faces by their country name
  using Face_ = Countries_arr::Face_handle::value_type;
  std::map<std::string, std::vector<Face_*>>  country_faces_map;
  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    auto& face = *fit;
    const auto& country_name = fit->data();
    // skipping spherical-face
    if (country_name.empty())
      continue;
    country_faces_map[country_name].push_back(&face);
  }

  std::cout << "triangulating individual faces\n";
  Country_triangles_map result;
  for (auto& [country_name, faces] : country_faces_map) {
    // CONVERT the face-points to QVector3D
    using Face_points = std::vector<QVector3D>;
    using Faces_ = std::vector<Face_points>;
    Faces_  all_faces_of_current_country;
    for (auto* face : faces) {
      // skip any face with no OUTER-CCB
      if (0 == face->number_of_outer_ccbs()) continue;

      std::vector<QVector3D> face_points;
      // loop on the egdes of the current outer-ccb
      auto first = face->outer_ccb();
      auto curr = first;
      do {
        auto ap = approx(curr->source()->point());
        QVector3D p(ap.dx(), ap.dy(), ap.dz());
        p.normalize();
        face_points.push_back(p);
      } while (++curr != first);

      all_faces_of_current_country.push_back(std::move(face_points));
    }

    // RESULTING TRIANGLE POINTS (every 3 point => triangle)
    auto& triangles = result[country_name];
    // loop on all approximated faces
    for (auto& face_points : all_faces_of_current_country) {
      //std::cout << "num face points = " << face_points.size() << std::endl;
      // no need to triangulate if the number of points is 3
      if (face_points.size() == 3) {
        triangles.insert(triangles.end(), face_points.begin(),
                         face_points.end());
          continue;
      }

      // find the centroid of all face-points
      QVector3D centroid(0, 0, 0);
      for (const auto& fp : face_points)
        centroid += fp;
      centroid /= face_points.size();
      centroid.normalize();
      auto normal = centroid;

      K::Point_3  plane_origin(centroid.x(), centroid.y(), centroid.z());
      K::Vector_3 plane_normal(normal.x(), normal.y(), normal.z());
      K::Plane_3 plane(plane_origin, plane_normal);

      Polygon_2 polygon;

      // project all points onto the plane
      K::Point_3 origin(0, 0, 0);
      for (const auto& fp : face_points) {
        // define a ray through the origin and the current point
        K::Point_3 current_point(fp.x(), fp.y(), fp.z());
        K::Ray_3 ray(origin, current_point);

        auto intersection = CGAL::intersection(plane, ray);
        if (!intersection.has_value())
          std::cout << "INTERSECTION ASSERTION ERROR!!!\n";
        auto ip = std::get<K::Point_3>(intersection.value());
        auto ip2 = plane.to_2d(ip);

        // add this to the polygon constraint
        polygon.push_back(ip2);
      }

      CDT cdt;
      cdt.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(),
                            true);

      std::unordered_map<Face_handle, bool> in_domain_map;
      boost::associative_property_map< std::unordered_map<Face_handle, bool> >
        in_domain(in_domain_map);

      //Mark facets that are inside the domain bounded by the polygon
      CGAL::mark_domain_in_triangulation(cdt, in_domain);

      // loop on all the triangles ("faces" in triangulation doc)
      for (Face_handle f : cdt.finite_face_handles()) {
        // if the current triangles is not inside the polygon -> skip it
        if (false == get(in_domain, f)) continue;

        for (int i = 0; i < 3; ++i) {
          auto tp = f->vertex(i)->point();
          auto tp3 = plane.to_3d(tp);
          QVector3D p3(tp3.x(), tp3.y(), tp3.z());
          p3.normalize();
          triangles.push_back(p3);
        }
      }
    }
  }

  return result;
}
