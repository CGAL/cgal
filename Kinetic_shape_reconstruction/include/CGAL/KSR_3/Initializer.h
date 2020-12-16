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

#ifndef CGAL_KSR_3_INITIALIZER_H
#define CGAL_KSR_3_INITIALIZER_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR_3/Data_structure.h>
#include <CGAL/KSR_3/Polygon_splitter.h>

namespace CGAL {
namespace KSR_3 {

template<typename GeomTraits>
class Initializer {

public:
  using Kernel = GeomTraits;

private:
  using FT          = typename Kernel::FT;
  using Point_2     = typename Kernel::Point_2;
  using Point_3     = typename Kernel::Point_3;
  using Segment_3   = typename Kernel::Segment_3;
  using Transform_3 = typename Kernel::Aff_transformation_3;

  using Data_structure   = KSR_3::Data_structure<Kernel>;
  using Polygon_splitter = KSR_3::Polygon_splitter<Data_structure, Kernel>;

  using Bbox_3  = CGAL::Bbox_3;
  using IVertex = typename Data_structure::IVertex;

public:
  Initializer(
    const bool debug,
    const bool verbose) :
  m_debug(debug),
  m_verbose(verbose),
  m_data(m_debug)
  { }

  template<
  typename InputRange,
  typename PolygonMap>
  const double initialize(
    const InputRange& input_range,
    const PolygonMap polygon_map,
    const unsigned int k,
    const double enlarge_bbox_ratio,
    const bool reorient) {

    if (m_verbose) {
      std::cout << std::endl << "--- INITIALIZING KSR:" << std::endl;
    }

    FT time_step;
    std::array<Point_3, 8> bbox;
    create_bounding_box(
      input_range, polygon_map,
      static_cast<FT>(enlarge_bbox_ratio),
      reorient, bbox, time_step);
    if (m_verbose) {
      std::cout << "* precomputed time_step: " << time_step << std::endl;
    }

    std::vector< std::vector<Point_3> > bbox_faces;
    bounding_box_to_polygons(bbox, bbox_faces);
    add_polygons(input_range, polygon_map, bbox_faces);

    if (m_verbose) std::cout << "* intersecting input polygons ... ";
    if (m_debug) {
      KSR_3::dump(m_data, "init");
      // KSR_3::dump_segmented_edges(m_data, "init");
    }

    m_data.check_integrity();
    make_polygons_intersection_free();
    m_data.check_integrity();
    set_k_intersections(k);

    if (m_verbose) std::cout << "done" << std::endl;
    if (m_debug) {
      KSR_3::dump(m_data, "intersected");
      // KSR_3::dump_segmented_edges(m_data, "intersected");
    }

    // for (KSR::size_t i = 6; i < m_data.number_of_support_planes(); ++i) {
    //   const auto& sp = m_data.support_plane(i);
    //   std::cout << "plane index: " << i << std::endl;
    //   std::cout << "plane: " <<
    //   sp.plane().a() << ", " <<
    //   sp.plane().b() << ", " <<
    //   sp.plane().c() << ", " <<
    //   sp.plane().d() << std::endl;
    // }

    m_data.check_bbox();
    return CGAL::to_double(time_step);
  }

  template<typename DS>
  void convert(DS& ds) {

    ds.clear();
    m_data.convert(ds);
    m_data.clear();

    ds.check_integrity();
    ds.check_bbox();
  }

  void clear() {
    m_data.clear();
  }

private:
  const bool m_debug;
  const bool m_verbose;
  Data_structure m_data;

  template<
  typename InputRange,
  typename PolygonMap>
  void create_bounding_box(
    const InputRange& input_range,
    const PolygonMap polygon_map,
    const FT enlarge_bbox_ratio,
    const bool reorient,
    std::array<Point_3, 8>& bbox,
    FT& time_step) const {

    if (reorient)
      initialize_optimal_box(input_range, polygon_map, bbox);
    else
      initialize_axis_aligned_box(input_range, polygon_map, bbox);

    CGAL_assertion(bbox.size() == 8);
    time_step  = KSR::distance(bbox.front(), bbox.back());
    time_step /= FT(50);

    enlarge_bounding_box(enlarge_bbox_ratio, bbox);

    const auto& minp = bbox.front();
    const auto& maxp = bbox.back();
    if (m_verbose) {
      std::cout << "* bounding box minp: " <<
      minp.x() << "\t, " << minp.y() << "\t, " << minp.z() << std::endl;
    }
    if (m_verbose) {
      std::cout << "* bounding box maxp: " <<
      maxp.x() << "\t, " << maxp.y() << "\t, " << maxp.z() << std::endl;
    }
  }

  template<
  typename InputRange,
  typename PolygonMap>
  void initialize_optimal_box(
    const InputRange& input_range,
    const PolygonMap polygon_map,
    std::array<Point_3, 8>& bbox) const {

    CGAL_assertion_msg(false, "TODO: IMPLEMENT THE ORIENTED OPTIMAL BBOX!");
  }

  template<
  typename InputRange,
  typename PolygonMap>
  void initialize_axis_aligned_box(
    const InputRange& input_range,
    const PolygonMap polygon_map,
    std::array<Point_3, 8>& bbox) const {

    Bbox_3 box;
    for (const auto& item : input_range) {
      const auto& polygon = get(polygon_map, item);
      box += CGAL::bbox_3(polygon.begin(), polygon.end());
    }

    // The order of faces corresponds to the standard order from here:
    // https://doc.cgal.org/latest/BGL/group__PkgBGLHelperFct.html#gad9df350e98780f0c213046d8a257358e
    bbox = {
      Point_3(box.xmin(), box.ymin(), box.zmin()),
      Point_3(box.xmax(), box.ymin(), box.zmin()),
      Point_3(box.xmax(), box.ymax(), box.zmin()),
      Point_3(box.xmin(), box.ymax(), box.zmin()),
      Point_3(box.xmin(), box.ymax(), box.zmax()),
      Point_3(box.xmin(), box.ymin(), box.zmax()),
      Point_3(box.xmax(), box.ymin(), box.zmax()),
      Point_3(box.xmax(), box.ymax(), box.zmax()) };
  }

  void enlarge_bounding_box(
    const FT enlarge_bbox_ratio,
    std::array<Point_3, 8>& bbox) const {

    CGAL_assertion_msg(
      enlarge_bbox_ratio > FT(1), "TODO: HANDLE THE CASE ENLARGE_BBOX_RATIO = FT(1)");
    const auto a = CGAL::centroid(bbox.begin(), bbox.end());
    Transform_3 scale(CGAL::Scaling(), enlarge_bbox_ratio);
    for (auto& point : bbox)
      point = scale.transform(point);

    const auto b = CGAL::centroid(bbox.begin(), bbox.end());
    Transform_3 translate(CGAL::Translation(), a - b);
    for (auto& point : bbox)
      point = translate.transform(point);
  }

  void bounding_box_to_polygons(
    const std::array<Point_3, 8>& bbox,
    std::vector< std::vector<Point_3> >& bbox_faces) const {

    bbox_faces.clear();
    bbox_faces.reserve(6);

    bbox_faces.push_back({bbox[0], bbox[1], bbox[2], bbox[3]});
    bbox_faces.push_back({bbox[0], bbox[1], bbox[6], bbox[5]});
    bbox_faces.push_back({bbox[1], bbox[2], bbox[7], bbox[6]});
    bbox_faces.push_back({bbox[2], bbox[3], bbox[4], bbox[7]});
    bbox_faces.push_back({bbox[3], bbox[0], bbox[5], bbox[4]});
    bbox_faces.push_back({bbox[5], bbox[6], bbox[7], bbox[4]});
    CGAL_assertion(bbox_faces.size() == 6);

    // Simon's bbox. The faces are different.
    // const FT xmin = bbox[0].x();
    // const FT ymin = bbox[0].y();
    // const FT zmin = bbox[0].z();
    // const FT xmax = bbox[7].x();
    // const FT ymax = bbox[7].y();
    // const FT zmax = bbox[7].z();
    // const std::vector<Point_3> sbbox = {
    //   Point_3(xmin, ymin, zmin),
    //   Point_3(xmin, ymin, zmax),
    //   Point_3(xmin, ymax, zmin),
    //   Point_3(xmin, ymax, zmax),
    //   Point_3(xmax, ymin, zmin),
    //   Point_3(xmax, ymin, zmax),
    //   Point_3(xmax, ymax, zmin),
    //   Point_3(xmax, ymax, zmax) };

    // bbox_faces.push_back({sbbox[0], sbbox[1], sbbox[3], sbbox[2]});
    // bbox_faces.push_back({sbbox[4], sbbox[5], sbbox[7], sbbox[6]});
    // bbox_faces.push_back({sbbox[0], sbbox[1], sbbox[5], sbbox[4]});
    // bbox_faces.push_back({sbbox[2], sbbox[3], sbbox[7], sbbox[6]});
    // bbox_faces.push_back({sbbox[1], sbbox[5], sbbox[7], sbbox[3]});
    // bbox_faces.push_back({sbbox[0], sbbox[4], sbbox[6], sbbox[2]});
    // CGAL_assertion(bbox_faces.size() == 6);
  }

  template<
  typename InputRange,
  typename PolygonMap>
  void add_polygons(
    const InputRange& input_range,
    const PolygonMap polygon_map,
    const std::vector< std::vector<Point_3> >& bbox_faces) {

    m_data.reserve(input_range.size());
    add_bbox_faces(bbox_faces);
    add_input_polygons(input_range, polygon_map);
  }

  void add_bbox_faces(
    const std::vector< std::vector<Point_3> >& bbox_faces) {

    for (const auto& bbox_face : bbox_faces)
      m_data.add_bbox_polygon(bbox_face);

    CGAL_assertion(m_data.number_of_support_planes() == 6);
    CGAL_assertion(m_data.ivertices().size() == 8);
    CGAL_assertion(m_data.iedges().size() == 12);

    if (m_verbose) {
      std::cout << "* inserted bbox faces: " << bbox_faces.size() << std::endl;
    }
  }

  template<
  typename InputRange,
  typename PolygonMap>
  void add_input_polygons(
    const InputRange& input_range,
    const PolygonMap polygon_map) {

    KSR::size_t input_index = 0;
    for (const auto& item : input_range) {
      const auto& polygon = get(polygon_map, item);
      m_data.add_input_polygon(polygon, input_index);
      ++input_index;
    }
    CGAL_assertion(m_data.number_of_support_planes() > 6);
    if (m_verbose) {
      std::cout << "* inserted input polygons: " <<  input_range.size() << std::endl;
    }
  }

  void make_polygons_intersection_free() {

    // First, create all transverse intersection lines.
    using Map_p2vv = std::map<KSR::Idx_set, std::pair<IVertex, IVertex> >;
    Map_p2vv map_p2vv;

    for (const auto ivertex : m_data.ivertices()) {
      const auto key = m_data.intersected_planes(ivertex, false);
      if (key.size() < 2) {
        continue;
      }

      const auto pair = map_p2vv.insert(
        std::make_pair(key, std::make_pair(ivertex, IVertex())));
      const bool is_inserted = pair.second;
      if (!is_inserted) {
        pair.first->second.second = ivertex;
      }
    }

    // Then, intersect these lines to find internal intersection vertices.
    using Pair_pv = std::pair< KSR::Idx_set, KSR::vector<IVertex> >;
    KSR::vector<Pair_pv> todo;
    for (auto it_a = map_p2vv.begin(); it_a != map_p2vv.end(); ++it_a) {
      const auto& set_a = it_a->first;

      todo.push_back(std::make_pair(set_a, KSR::vector<IVertex>()));
      auto& crossed_vertices = todo.back().second;
      crossed_vertices.push_back(it_a->second.first);

      std::set<KSR::Idx_set> done;
      for (auto it_b = map_p2vv.begin(); it_b != map_p2vv.end(); ++it_b) {
        const auto& set_b = it_b->first;

        KSR::size_t common_plane_idx = KSR::no_element();
        std::set_intersection(
          set_a.begin(), set_a.end(), set_b.begin(), set_b.end(),
          boost::make_function_output_iterator(
            [&](const KSR::size_t idx) -> void {
              common_plane_idx = idx;
            }
          )
        );

        if (common_plane_idx != KSR::no_element()) {
          auto union_set = set_a;
          union_set.insert(set_b.begin(), set_b.end());
          if (!done.insert(union_set).second) {
            continue;
          }

          Point_2 inter;
          if (!KSR::intersection(
            m_data.to_2d(common_plane_idx,
              Segment_3(m_data.point_3(it_a->second.first), m_data.point_3(it_a->second.second))),
            m_data.to_2d(common_plane_idx,
              Segment_3(m_data.point_3(it_b->second.first), m_data.point_3(it_b->second.second))),
            inter)) {

            continue;
          }

          crossed_vertices.push_back(
            m_data.add_ivertex(m_data.to_3d(common_plane_idx, inter), union_set));
        }
      }
      crossed_vertices.push_back(it_a->second.second);
    }

    for (auto& t : todo) {
      m_data.add_iedge(t.first, t.second);
    }

    // Refine polygons.
    for (KSR::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      Polygon_splitter splitter(m_data);
      splitter.split_support_plane(i);
      // if (i >= 6 && m_debug) {
      //   KSR_3::dump(m_data, "intersected-iter-" + std::to_string(i));
      // }
    }
  }

  void set_k_intersections(const unsigned int k) {

    for (KSR::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      for (const auto pface : m_data.pfaces(i)) {
        m_data.k(pface) = k;
      }
    }
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_INITIALIZER_H
