// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_INITIALIZER_H
#define CGAL_KSR_3_INITIALIZER_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR_3/Data_structure.h>
#include <CGAL/KSR_3/Polygon_splitter.h>

namespace CGAL {
namespace KSR_3 {

// TODO: DOES NOT WORK WITH INEXACT KERNEL!
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

  using IVertex  = typename Data_structure::IVertex;
  using IK       = CGAL::Exact_predicates_inexact_constructions_kernel;
  using IFT      = typename IK::FT;
  using IPoint_3 = typename IK::Point_3;

  using Bbox_3     = CGAL::Bbox_3;
  using OBB_traits = CGAL::Oriented_bounding_box_traits_3<IK>;

public:
  Initializer(
    const bool verbose, const bool dprint, const bool debug, Data_structure& data) :
  m_verbose(verbose), m_export(dprint), m_debug(debug), m_data(data)
  { }

  template<
  typename InputRange,
  typename PolygonMap>
  double initialize(
    const InputRange& input_range,
    const PolygonMap polygon_map,
    const unsigned int k,
    const double enlarge_bbox_ratio,
    const bool reorient) {

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

    CGAL_assertion(m_data.check_integrity(false));
    make_polygons_intersection_free();
    CGAL_assertion(m_data.check_integrity(false));
    set_k_intersections(k);

    if (m_verbose) std::cout << "done" << std::endl;
    if (m_debug) {
      KSR_3::dump(m_data, "intersected");
      // KSR_3::dump_segmented_edges(m_data, "intersected");
    }

    // for (std::size_t i = 6; i < m_data.number_of_support_planes(); ++i) {
    //   const auto& sp = m_data.support_plane(i);
    //   std::cout << "plane index: " << i << std::endl;
    //   std::cout << "plane: " <<
    //   sp.plane().a() << ", " <<
    //   sp.plane().b() << ", " <<
    //   sp.plane().c() << ", " <<
    //   sp.plane().d() << std::endl;
    // }

    CGAL_assertion(m_data.check_bbox());
    m_data.set_limit_lines();
    m_data.precompute_iedge_data();
    CGAL_assertion(m_data.check_integrity());
    CGAL_assertion(m_data.check_intersection_graph());

    return CGAL::to_double(time_step);
  }

  template<typename DS>
  void transfer_to(DS& ds) {

    CGAL_assertion_msg(false,
    "USE THIS ONLY FOR CONVERTING EXACT DATA TO INEXACT DS!");
    ds.clear();
    m_data.transfer_to(ds);
    m_data.clear();

    CGAL_assertion(ds.check_integrity(false));
    CGAL_assertion(ds.check_bbox());
  }

  void clear() {
    // to be added
  }

private:
  const bool m_verbose;
  const bool m_export;
  const bool m_debug;
  Data_structure& m_data;

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

    if (reorient) {
      initialize_optimal_box(input_range, polygon_map, bbox);
    } else {
      initialize_axis_aligned_box(input_range, polygon_map, bbox);
    }

    CGAL_assertion(bbox.size() == 8);
    time_step  = KSR::distance(bbox.front(), bbox.back());
    time_step /= FT(50);

    enlarge_bounding_box(enlarge_bbox_ratio, bbox);

    const auto& minp = bbox.front();
    const auto& maxp = bbox.back();
    if (m_verbose) {
      std::cout.precision(20);
      std::cout << "* bounding box minp: " << std::fixed <<
      minp.x() << "\t, " << minp.y() << "\t, " << minp.z() << std::endl;
    }
    if (m_verbose) {
      std::cout.precision(20);
      std::cout << "* bounding box maxp: " << std::fixed <<
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

    // Number of input points.
    std::size_t num_points = 0;
    for (const auto& item : input_range) {
      const auto& polygon = get(polygon_map, item);
      num_points += polygon.size();
    }

    // Set points.
    std::vector<IPoint_3> ipoints;
    ipoints.reserve(num_points);
    for (const auto& item : input_range) {
      const auto& polygon = get(polygon_map, item);
      for (const auto& point : polygon) {
        const IPoint_3 ipoint(
          static_cast<IFT>(CGAL::to_double(point.x())),
          static_cast<IFT>(CGAL::to_double(point.y())),
          static_cast<IFT>(CGAL::to_double(point.z())));
        ipoints.push_back(ipoint);
      }
    }

    // Compute optimal bbox.
    // The order of faces corresponds to the standard order from here:
    // https://doc.cgal.org/latest/BGL/group__PkgBGLHelperFct.html#gad9df350e98780f0c213046d8a257358e
    const OBB_traits obb_traits;
    std::array<IPoint_3, 8> ibbox;
    CGAL::oriented_bounding_box(
      ipoints, ibbox,
      CGAL::parameters::use_convex_hull(true).
      geom_traits(obb_traits));

    for (std::size_t i = 0; i < 8; ++i) {
      const auto& ipoint = ibbox[i];
      const Point_3 point(
        static_cast<FT>(ipoint.x()),
        static_cast<FT>(ipoint.y()),
        static_cast<FT>(ipoint.z()));
      bbox[i] = point;
    }

    const FT bbox_length_1 = KSR::distance(bbox[0], bbox[1]);
    const FT bbox_length_2 = KSR::distance(bbox[0], bbox[3]);
    const FT bbox_length_3 = KSR::distance(bbox[0], bbox[5]);
    CGAL_assertion(bbox_length_1 >= FT(0));
    CGAL_assertion(bbox_length_2 >= FT(0));
    CGAL_assertion(bbox_length_3 >= FT(0));
    const FT tol = KSR::tolerance<FT>();
    if (bbox_length_1 < tol || bbox_length_2 < tol || bbox_length_3 < tol) {
      if (m_verbose) {
        std::cout << "* warning: optimal bounding box is flat, reverting ..." << std::endl;
      }
      initialize_axis_aligned_box(input_range, polygon_map, bbox);
    } else {
      if (m_verbose) {
        std::cout << "* using optimal bounding box" << std::endl;
      }
    }
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

    const FT bbox_length_1 = KSR::distance(bbox[0], bbox[1]);
    const FT bbox_length_2 = KSR::distance(bbox[0], bbox[3]);
    const FT bbox_length_3 = KSR::distance(bbox[0], bbox[5]);
    CGAL_assertion(bbox_length_1 >= FT(0));
    CGAL_assertion(bbox_length_2 >= FT(0));
    CGAL_assertion(bbox_length_3 >= FT(0));
    const FT tol = KSR::tolerance<FT>();
    if (bbox_length_1 < tol || bbox_length_2 < tol || bbox_length_3 < tol) {
      if (bbox_length_1 < tol) {

        CGAL_assertion_msg(bbox_length_2 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");
        CGAL_assertion_msg(bbox_length_3 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");
        const FT x = FT(40) * tol; // 40 is a magic number but it is not a big deal in this case

        bbox[0] = Point_3(bbox[0].x() - x, bbox[0].y(), bbox[0].z());
        bbox[3] = Point_3(bbox[3].x() - x, bbox[3].y(), bbox[3].z());
        bbox[4] = Point_3(bbox[4].x() - x, bbox[4].y(), bbox[4].z());
        bbox[5] = Point_3(bbox[5].x() - x, bbox[5].y(), bbox[5].z());

        bbox[1] = Point_3(bbox[1].x() + x, bbox[1].y(), bbox[1].z());
        bbox[2] = Point_3(bbox[2].x() + x, bbox[2].y(), bbox[2].z());
        bbox[7] = Point_3(bbox[7].x() + x, bbox[7].y(), bbox[7].z());
        bbox[6] = Point_3(bbox[6].x() + x, bbox[6].y(), bbox[6].z());
        if (m_verbose) std::cout << "* setting x-based flat axis-aligned bbox" << std::endl;

      } else if (bbox_length_2 < tol) {

        CGAL_assertion_msg(bbox_length_3 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");
        CGAL_assertion_msg(bbox_length_1 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");
        const FT y = FT(40) * tol; // 40 is a magic number but it is not a big deal in this case

        bbox[0] = Point_3(bbox[0].x(), bbox[0].y() - y, bbox[0].z());
        bbox[1] = Point_3(bbox[1].x(), bbox[1].y() - y, bbox[1].z());
        bbox[6] = Point_3(bbox[6].x(), bbox[6].y() - y, bbox[6].z());
        bbox[5] = Point_3(bbox[5].x(), bbox[5].y() - y, bbox[5].z());

        bbox[3] = Point_3(bbox[3].x(), bbox[3].y() + y, bbox[3].z());
        bbox[2] = Point_3(bbox[2].x(), bbox[2].y() + y, bbox[2].z());
        bbox[7] = Point_3(bbox[7].x(), bbox[7].y() + y, bbox[7].z());
        bbox[4] = Point_3(bbox[4].x(), bbox[4].y() + y, bbox[4].z());
        if (m_verbose) std::cout << "* setting y-based flat axis-aligned bbox" << std::endl;

      } else if (bbox_length_3 < tol) {

        CGAL_assertion_msg(bbox_length_1 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");
        CGAL_assertion_msg(bbox_length_2 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");
        const FT z = FT(40) * tol; // 40 is a magic number but it is not a big deal in this case

        bbox[0] = Point_3(bbox[0].x(), bbox[0].y(), bbox[0].z() - z);
        bbox[1] = Point_3(bbox[1].x(), bbox[1].y(), bbox[1].z() - z);
        bbox[2] = Point_3(bbox[2].x(), bbox[2].y(), bbox[2].z() - z);
        bbox[3] = Point_3(bbox[3].x(), bbox[3].y(), bbox[3].z() - z);

        bbox[5] = Point_3(bbox[5].x(), bbox[5].y(), bbox[5].z() + z);
        bbox[6] = Point_3(bbox[6].x(), bbox[6].y(), bbox[6].z() + z);
        bbox[7] = Point_3(bbox[7].x(), bbox[7].y(), bbox[7].z() + z);
        bbox[4] = Point_3(bbox[4].x(), bbox[4].y(), bbox[4].z() + z);
        if (m_verbose) std::cout << "* setting z-based flat axis-aligned bbox" << std::endl;

      } else {
        CGAL_assertion_msg(false, "ERROR: WRONG CASE!");
      }
    } else {
      if (m_verbose) {
        std::cout << "* using axis aligned bounding box" << std::endl;
      }
    }
  }

  void enlarge_bounding_box(
    const FT enlarge_bbox_ratio,
    std::array<Point_3, 8>& bbox) const {

    FT enlarge_ratio = enlarge_bbox_ratio;
    const FT tol = KSR::tolerance<FT>();
    if (enlarge_bbox_ratio == FT(1)) {
      enlarge_ratio += FT(2) * tol;
    }

    const auto a = CGAL::centroid(bbox.begin(), bbox.end());
    Transform_3 scale(CGAL::Scaling(), enlarge_ratio);
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

    using Polygon_2 = std::vector<Point_2>;
    using Indices = std::vector<std::size_t>;

    std::map< std::size_t, std::pair<Polygon_2, Indices> > polygons;
    preprocess_polygons(input_range, polygon_map, polygons);
    CGAL_assertion(polygons.size() > 0);

    for (const auto& item : polygons) {
      const std::size_t support_plane_idx = item.first;
      const auto& pair = item.second;
      const Polygon_2& polygon = pair.first;
      const Indices& input_indices = pair.second;
      m_data.add_input_polygon(support_plane_idx, input_indices, polygon);
    }

    CGAL_assertion(m_data.number_of_support_planes() > 6);
    if (m_verbose) {
      std::cout << "* inserted input polygons: " << polygons.size() << std::endl;
    }
  }

  template<typename PointRange>
  void convert_polygon(
    const std::size_t support_plane_idx,
    const PointRange& polygon_3,
    std::vector<Point_2>& polygon_2) {

    polygon_2.clear();
    polygon_2.reserve(polygon_3.size());
    for (const auto& point : polygon_3) {
      const Point_3 converted(
        static_cast<FT>(point.x()),
        static_cast<FT>(point.y()),
        static_cast<FT>(point.z()));
      polygon_2.push_back(
        m_data.support_plane(support_plane_idx).to_2d(converted));
    }
    CGAL_assertion(polygon_2.size() == polygon_3.size());
  }

  template<
  typename InputRange,
  typename PolygonMap>
  void preprocess_polygons(
    const InputRange& input_range,
    const PolygonMap polygon_map,
    std::map< std::size_t, std::pair<
      std::vector<Point_2>,
      std::vector<std::size_t> > >& polygons) {

    std::size_t input_index = 0;
    std::vector<Point_2> polygon_2;
    std::vector<std::size_t> input_indices;
    for (const auto& item : input_range) {
      const auto& polygon_3 = get(polygon_map, item);

      bool is_added = true;
      std::size_t support_plane_idx = KSR::no_element();
      std::tie(support_plane_idx, is_added) = m_data.add_support_plane(polygon_3);
      CGAL_assertion(support_plane_idx != KSR::no_element());
      convert_polygon(support_plane_idx, polygon_3, polygon_2);

      if (is_added) {
        input_indices.clear();
        input_indices.push_back(input_index);
        polygons[support_plane_idx] = std::make_pair(polygon_2, input_indices);
      } else {
        CGAL_assertion_msg(false, "TODO: FINISH POLYGONS PREPROCESSING!");
      }
      ++input_index;
    }
  }

  void make_polygons_intersection_free() {

    if (m_debug) {
      std::cout << std::endl;
      std::cout.precision(20);
    }

    // First, create all transverse intersection lines.
    using Map_p2vv = std::map<std::set<std::size_t>, std::pair<IVertex, IVertex> >;
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
    using Pair_pv = std::pair< std::set<std::size_t>, std::vector<IVertex> >;
    std::vector<Pair_pv> todo;
    for (auto it_a = map_p2vv.begin(); it_a != map_p2vv.end(); ++it_a) {
      const auto& set_a = it_a->first;

      todo.push_back(std::make_pair(set_a, std::vector<IVertex>()));
      auto& crossed_vertices = todo.back().second;
      crossed_vertices.push_back(it_a->second.first);

      std::set<std::set<std::size_t>> done;
      for (auto it_b = map_p2vv.begin(); it_b != map_p2vv.end(); ++it_b) {
        const auto& set_b = it_b->first;

        std::size_t common_plane_idx = KSR::no_element();
        std::set_intersection(
          set_a.begin(), set_a.end(), set_b.begin(), set_b.end(),
          boost::make_function_output_iterator(
            [&](const std::size_t idx) -> void {
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

          Point_2 point;
          if (!KSR::intersection(
            m_data.to_2d(common_plane_idx,
              Segment_3(m_data.point_3(it_a->second.first), m_data.point_3(it_a->second.second))),
            m_data.to_2d(common_plane_idx,
              Segment_3(m_data.point_3(it_b->second.first), m_data.point_3(it_b->second.second))),
            point)) {

            continue;
          }

          crossed_vertices.push_back(
            m_data.add_ivertex(m_data.to_3d(common_plane_idx, point), union_set));
        }
      }
      crossed_vertices.push_back(it_a->second.second);
    }

    for (auto& t : todo) {
      m_data.add_iedge(t.first, t.second);
    }

    // Refine polygons.
    for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      Polygon_splitter splitter(m_data);
      splitter.split_support_plane(i);
      // if (i >= 6 && m_export) {
      //   KSR_3::dump(m_data, "intersected-iter-" + std::to_string(i));
      // }
    }
    // exit(EXIT_SUCCESS);
  }

  void set_k_intersections(const unsigned int k) {

    for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      for (const auto pface : m_data.pfaces(i)) {
        m_data.k(pface) = k;
      }
    }
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_INITIALIZER_H
