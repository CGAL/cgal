// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
// #include <CGAL/draw_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Arr_batched_point_location.h>

#include "Aos_defs.h"
#include "Aos_triangulator.h"

using json = nlohmann::ordered_json;

namespace {

  // use this traits every time you construct an arrangement!
  static Geom_traits s_traits;
  using Dir3 = Kernel::Direction_3;
  using Approximate_number_type = Geom_traits::Approximate_number_type;
  using Approximate_kernel = Geom_traits::Approximate_kernel;
  using Approximate_Vector_3 = CGAL::Vector_3<Approximate_kernel>;
  using Approximate_Direction_3 = Approximate_kernel::Direction_3;
  using Direction_3 = Kernel::Direction_3;
}

//! \brief
std::vector<QVector3D> Aos_triangulator::get_all(Aos::Arr_handle arrh) {
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Vb = CGAL::Triangulation_vertex_base_2<K>;
  using Fb = CGAL::Constrained_triangulation_face_base_2<K>;
  using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Itag = CGAL::Exact_predicates_tag;
  using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>;
  using Face_handle = CDT::Face_handle;
  using Polygon_2 = CGAL::Polygon_2<K>;

  auto& arr = *reinterpret_cast<Arrangement*>(arrh.get());

  //Geom_traits traits;
  auto approx = s_traits.approximate_2_object();

  std::vector<std::vector<QVector3D>> all_faces;
  // loop on all faces of the arrangement
  for (auto fh : arr.face_handles()) {
    // skip any face with no OUTER-CCB
    if (0 == fh->number_of_outer_ccbs()) continue;

    // COMPUTE THE CENTROID OF ALL FACE-POINTS
    std::vector<QVector3D> face_points;

    // loop on the edges of the current outer-ccb
    auto first = fh->outer_ccb();
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
  std::vector<QVector3D> triangles;

  std::cout << "triangulating individual faces\n";

  // loop on all approximated faces
  for (auto& face_points : all_faces) {
    std::cout << "num face points = " << face_points.size() << std::endl;
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
      if (! get(in_domain, f)) continue;

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
Country_triangles_map
Aos_triangulator::get_by_country(Aos::Arr_handle arrh, float error,
                                 std::size_t num_uniform_points) {

  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Pnt_2 = typename K::Point_2;
  using Pnt_3 = typename K::Point_3;
  using Vb = CGAL::Triangulation_vertex_base_2<K>;
  using Fb = CGAL::Constrained_triangulation_face_base_2<K>;
  using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Itag = CGAL::Exact_predicates_tag;
  using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>;
  using Face_handle = CDT::Face_handle;
  using Approximate_point_2 = Geom_traits::Approximate_point_2;
  using Aos = Countries_arr;
  using Aos_pnt = Aos::Point_2;

  auto& arr = *reinterpret_cast<Aos*>(arrh.get());
  auto approx = s_traits.approximate_2_object();

  // group the faces by their country name
  using Face_const_handle = Countries_arr::Face_const_handle;
  std::map<std::string, std::vector<Face_const_handle>> country_faces_map;
  for (auto fh : arr.face_handles()) {
    const auto& country_name = fh->data();
    // skipping spherical-face
    if (country_name.empty()) continue;
    country_faces_map[country_name].push_back(fh);
  }

  // Generate random point uniformly distributed
  std::vector<Pnt_3> rps;
  rps.reserve(num_uniform_points);
  CGAL::Random_points_in_sphere_3<Pnt_3> g(1.0);
  std::copy_n(g, rps.capacity(), std::back_inserter(rps));
  std::map<Face_const_handle, std::vector<Pnt_3>> face_random_points;
  using Point_location_result = CGAL::Arr_point_location_result<Aos>;
  using Query_result = std::pair<Aos_pnt, Point_location_result::Type>;
  std::vector<Aos_pnt> rqs(rps.size());
  auto ctr_p = arr.geometry_traits()->construct_point_2_object();
  std::transform(rps.begin(), rps.end(), rqs.begin(),
                 [&](const Pnt_3& rp) -> Aos_pnt {
                   return ctr_p(rp.x(), rp.y(), rp.z());
                 });
  std::list<Query_result> results;
  CGAL::locate(arr, rqs.begin(), rqs.end(), std::back_inserter(results));
  for (auto& [q, res] : results) {
    const Aos::Face_const_handle* f;
    if ((f = std::get_if<Face_const_handle>(&res))) {
      auto x = CGAL::to_double(q.dx());
      auto y = CGAL::to_double(q.dy());
      auto z = CGAL::to_double(q.dz());
      face_random_points[*f].push_back(Pnt_3(x, y, z));
    }
  }

  // Triangulate the faces
  Country_triangles_map result;
  for (auto& [country_name, fhs] : country_faces_map) {
    // std::cout << "processing country " << country_name << std::endl;
    auto& triangles = result[country_name];
    // CONVERT the face-points to QVector3D
    for (auto fh : fhs) {
      // skip any face with no OUTER-CCB
      if (0 == fh->number_of_outer_ccbs()) continue;

      std::vector<QVector3D> face_points;
      // Loop on the edges of the current outer-ccb
      auto first = fh->outer_ccb();
      auto curr = first;
      do {
        std::vector<Approximate_point_2> apx_points;
        const auto& xcv = curr->curve();
        approx(xcv, error, std::back_insert_iterator(apx_points),
               curr->direction() == CGAL::ARR_LEFT_TO_RIGHT);
        for (auto it = apx_points.begin(); it != apx_points.end()-1; ++it) {
          const auto& apx_p = *it;
          const QVector3D p(apx_p.dx(), apx_p.dy(), apx_p.dz());
          face_points.push_back(p);
        }
      } while (++curr != first);

      // no need to triangulate if the number of points is 3
      if (face_points.size() == 3) {
        triangles.insert(triangles.end(), face_points.begin(),
                         face_points.end());
          continue;
      }

      // Calculate the centroid of all face-points
      QVector3D centroid(0, 0, 0);
      for (const auto& fp : face_points) centroid += fp;
      centroid /= face_points.size();
      centroid.normalize();
      auto normal = centroid;

      Pnt_3 plane_origin(centroid.x(), centroid.y(), centroid.z());
      K::Vector_3 plane_normal(normal.x(), normal.y(), normal.z());
      K::Plane_3 plane(plane_origin, plane_normal);

      // Project all points onto the plane
      std::vector<Pnt_2> planar_points(face_points.size());
      std::transform(face_points.begin(), face_points.end(),
                     planar_points.begin(),
                     [&](const QVector3D& fp) -> Pnt_2 {
                       // define a ray through the origin and the current point
                       Pnt_3 p3(fp.x(), fp.y(), fp.z());
                       K::Ray_3 ray(CGAL::ORIGIN, p3);
                       auto intersection = CGAL::intersection(plane, ray);
                       if (! intersection.has_value())
                         std::cout << "INTERSECTION ASSERTION ERROR!!!\n";
                       auto ip = std::get<Pnt_3>(intersection.value());
                       return plane.to_2d(ip);
                     });
      CDT cdt;

      // Insert points uniformly distributed into the triangulation
      auto it = face_random_points.find(fh);
      if (it != face_random_points.end()) {
        const auto& points = it->second;
        for (const auto& p3 : points) {
          K::Ray_3 ray(CGAL::ORIGIN, p3);
          auto intersection = CGAL::intersection(plane, ray);
          if (! intersection.has_value())
            std::cout << "INTERSECTION ASSERTION ERROR!!!\n";
          auto ip = std::get<Pnt_3>(intersection.value());
          auto p2 = plane.to_2d(ip);
          cdt.insert(p2);
        }
      }

      // Insert the constraints into the triangulation
      cdt.insert_constraint(planar_points.begin(), planar_points.end(), true);

      // Mark facets that are inside the domain bounded by the polygon
      std::unordered_map<Face_handle, bool> in_domain_map;
      boost::associative_property_map< std::unordered_map<Face_handle, bool>>
        in_domain(in_domain_map);
      CGAL::mark_domain_in_triangulation(cdt, in_domain);

      // Loop on all the triangles ("faces" in triangulation doc)
      for (Face_handle f : cdt.finite_face_handles()) {
        // If the current triangles is not inside the polygon -> skip it
        if (! get(in_domain, f)) continue;

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
