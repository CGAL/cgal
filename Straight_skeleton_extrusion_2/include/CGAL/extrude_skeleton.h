// Copyright (c) 2023 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//

#ifndef CGAL_SLS_EXTRUSION_EXTRUDE_SKELETON_H
#define CGAL_SLS_EXTRUSION_EXTRUDE_SKELETON_H

#include <CGAL/license/Straight_skeleton_extrusion_2.h>

#ifdef CGAL_SLS_DEBUG_DRAW
  #include <CGAL/Straight_skeleton_2/IO/print.h>
  #include <CGAL/draw_straight_skeleton_2.h>
  #include <CGAL/draw_polygon_2.h>
  #include <CGAL/draw_polygon_with_holes_2.h>
  #include <CGAL/draw_triangulation_2.h>
#endif

#include <CGAL/create_weighted_straight_skeleton_2.h>
#include <CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/create_weighted_offset_polygons_2.h>
#include <CGAL/create_weighted_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_converter_2.h>

#include <CGAL/General_polygon_with_holes_2.h>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#ifdef CGAL_SLS_OUTPUT_FILES
#include <CGAL/IO/polygon_soup_io.h>
#endif

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/mark_domain_in_triangulation.h>

#include <CGAL/boost/graph/border.h>
#include <CGAL/Cartesian/function_objects.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/enum.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/unordered_flat_map.h>
#include <CGAL/utility.h>

#include <boost/algorithm/clamp.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/range/value_type.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace CGAL {
namespace Straight_skeleton_extrusion {
namespace internal {

template <typename K>
static constexpr bool is_exact_kernel =
  std::is_same_v<typename CGAL::Algebraic_structure_traits<typename K::FT>::Is_exact, CGAL::Tag_true>;

template <typename SK>
struct Extrusion_kernel_selector
{
  using type = std::conditional_t<is_exact_kernel<SK>, SK, CGAL::Exact_predicates_exact_constructions_kernel>;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Below is to handle vertical slabs.

enum class Slope
{
  UNKNOWN = 0,
  INWARD,
  OUTWARD,
  VERTICAL
};

// Vertical slabs cannot be represented with (strictly) positive weights.
// The workaround is to give a very large weight (see preprocess_weights()), which creates
// near vertical faces.
// If the macro below is defined, the points of these faces are moved onto the vertical planes
// after the output mesh has been constructed.
#define CGAL_SLS_SNAP_TO_VERTICAL_SLABS

template <typename StraightSkeleton_2,
          typename OffsetBuilderTraits,
          typename XK>
class Skeleton_offset_correspondence_builder_visitor
  : public CGAL::Default_polygon_offset_builder_2_visitor<OffsetBuilderTraits, StraightSkeleton_2>
{
  using Base = CGAL::Default_polygon_offset_builder_2_visitor<OffsetBuilderTraits, StraightSkeleton_2>;

  using FT = typename XK::FT;
  using Point_2 = typename XK::Point_2;
  using Point_3 = typename XK::Point_3;

  using SS_Halfedge_const_handle = typename StraightSkeleton_2::Halfedge_const_handle;

  using HDS = typename StraightSkeleton_2::Base;
  using HDS_Vertex_const_handle = typename HDS::Vertex_const_handle;
  using HDS_Halfedge_const_handle = typename HDS::Halfedge_const_handle;

public:
  Skeleton_offset_correspondence_builder_visitor(const FT& height,
                                                 std::vector<Point_3>& points,
                                                 CGAL::unordered_flat_map<HDS_Vertex_const_handle, std::size_t>& vertex_point_ids,
                                                 CGAL::unordered_flat_map<HDS_Halfedge_const_handle, std::size_t>& halfedge_point_ids,
                                                 std::vector<std::vector<HDS_Halfedge_const_handle> >& contour_halfedges,
                                                 const bool ignore_first_contour = false) // for outer skeletons
    : m_height(height),
      m_points(points),
      m_vertex_point_ids(vertex_point_ids),
      m_halfedge_point_ids(halfedge_point_ids),
      m_contour_halfedges(contour_halfedges),
      m_skip_next_contour(ignore_first_contour)
  { }

private:
  // Same as Polygon_offset_builder_2::Compare_offset_against_event_time(), plus contour nodes.
  CGAL::Comparison_result compare_height_to_node_time(HDS_Vertex_const_handle node) const
  {
    if(node->is_contour())
      return CGAL::LARGER; // abs height is > 0
    if(node->has_infinite_time())
      return CGAL::SMALLER;

    return static_cast<CGAL::Comparison_result>(
      OffsetBuilderTraits().compare_offset_against_event_time_2_object()(CGAL::abs(m_height),
                                                                         node->trisegment()));
  }

public:
  void on_offset_contour_started() const
  {
    // std::cout << "~~ new contour ~~" << std::endl;
    m_tentative_contour.clear();
  }

  void on_offset_point(const Point_2& op,
                       SS_Halfedge_const_handle hook) const
  {
    CGAL_precondition(hook->is_bisector());

    // std::cout << "-> on offset point(" << op << ", " << " on H" << hook->id() << ")" << std::endl;
    m_tentative_contour.emplace_back(op, hook);
  }

  void on_offset_contour_finished(const bool is_complete) const
  {
    // if we have to ignore the first contour, that means it's the outer frame's
    // and the outer frame should be the first contour, and a complete contour.
    CGAL_assertion(!m_skip_next_contour || is_complete);

    if(!is_complete)
      return;

    if(std::exchange(m_skip_next_contour, false))
      return;

    m_contour_halfedges.emplace_back();
    std::vector<HDS_Halfedge_const_handle>& contour = m_contour_halfedges.back();

     std::vector<std::size_t> ids; // as to detect consecutive hooks sharing the same point

     for(const auto& e : m_tentative_contour)
     {
        const Point_2& op = e.first;
        HDS_Halfedge_const_handle hook = e.second;
        HDS_Halfedge_const_handle canonical_hook = (hook < hook->opposite()) ? hook : hook->opposite();

        // the visitor could have given the hook position...
        HDS_Vertex_const_handle node;
        if(compare_height_to_node_time(hook->vertex()) == CGAL::EQUAL)
          node = hook->vertex();
        // else if, because the slope is non-zero for hooks
        else if(compare_height_to_node_time(hook->opposite()->vertex()) == CGAL::EQUAL)
          node = hook->opposite()->vertex();

        std::size_t id = m_points.size();
        if(node != HDS_Vertex_const_handle())
        {
          auto vres = m_vertex_point_ids.emplace(node, id);
          if(vres.second) // successful insertion, never investigated this vertex before
          {
            // can't be a known halfedge with hook at source/target and unknown vertex
            CGAL_assertion(m_halfedge_point_ids.count(canonical_hook) == 0);

            m_points.push_back(Point_3{op.x(), op.y(), m_height});
          }
          else
          {
            id = vres.first->second;
          }

          m_halfedge_point_ids.emplace(canonical_hook, id);
        }
        else
        {
          auto res = m_halfedge_point_ids.emplace(canonical_hook, id);
          if(res.second) // successful insertion, never investigated this halfedge before
            m_points.push_back(Point_3{op.x(), op.y(), m_height});
        }

        // Two consecutive hooks can share a node at the offset time (offset contouring removes such
        // duplicates geometrically, which is not reliable with inexact constructions)
        if(!ids.empty() && ids.back() == id)
          continue;

        ids.push_back(id);
        contour.push_back(canonical_hook);
     }

    if(ids.size() > 1 && ids.front() == ids.back())
      contour.pop_back();
   }

private:
  const FT& m_height;
  std::vector<Point_3>& m_points;
  CGAL::unordered_flat_map<HDS_Vertex_const_handle, std::size_t>& m_vertex_point_ids;
  CGAL::unordered_flat_map<HDS_Halfedge_const_handle, std::size_t>& m_halfedge_point_ids;
  std::vector<std::vector<HDS_Halfedge_const_handle> >& m_contour_halfedges;
  mutable std::vector<std::pair<Point_2, HDS_Halfedge_const_handle> > m_tentative_contour;
  mutable bool m_skip_next_contour;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename SK,
          typename XK>
class Extrusion_builder
{
  using SFT = typename SK::FT;
  using FT = typename XK::FT;
  using Point_2 = typename XK::Point_2;
  using Line_2 = typename XK::Line_2;
  using Point_3 = typename XK::Point_3;
  using Vector_3 = typename XK::Vector_3;

  using Polygon_2 = CGAL::Polygon_2<XK>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<XK>;

  using Offset_polygons = std::vector<std::shared_ptr<Polygon_2> >;
  using Offset_polygons_with_holes = std::vector<std::shared_ptr<Polygon_with_holes_2> >;

  using S_Straight_skeleton_2 = CGAL::Straight_skeleton_2<SK>;
  using S_Straight_skeleton_2_ptr = std::shared_ptr<S_Straight_skeleton_2>;
  using Straight_skeleton_2 = CGAL::Straight_skeleton_2<XK>;
  using Straight_skeleton_2_ptr = std::shared_ptr<Straight_skeleton_2>;

  using SS_Vertex_const_handle = typename Straight_skeleton_2::Vertex_const_handle;
  using SS_Halfedge_const_handle = typename Straight_skeleton_2::Halfedge_const_handle;

  using HDS = typename Straight_skeleton_2::Base;
  using HDS_Vertex_const_handle = typename HDS::Vertex_const_handle;
  using HDS_Halfedge_const_handle = typename HDS::Halfedge_const_handle;
  using HDS_Halfedge_const_iterator = typename HDS::Halfedge_const_iterator;
  using HDS_Face_handle = typename HDS::Face_handle;

  using Offset_builder_traits = CGAL::Polygon_offset_builder_traits_2<XK>;
  using Offset_builder_no_visitor = CGAL::Polygon_offset_builder_2<Straight_skeleton_2,
                                                                   Offset_builder_traits,
                                                                   Polygon_2>;
  using Visitor = Skeleton_offset_correspondence_builder_visitor<Straight_skeleton_2, Offset_builder_traits, XK>;
  using Offset_builder = CGAL::Polygon_offset_builder_2<Straight_skeleton_2,
                                                        Offset_builder_traits,
                                                        Polygon_2,
                                                        Visitor>;

private:
  SK m_sk;
  XK m_xk;

  bool m_verbose;

private:
  // these are reset before each run
  CGAL::unordered_flat_map<HDS_Vertex_const_handle, std::size_t> m_vertex_point_ids;
  CGAL::unordered_flat_map<HDS_Halfedge_const_handle, std::size_t> m_halfedge_point_ids;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
  // This is to deal with vertical slabs: we temporarily give a non-vertical slope to be able
  // to construct the skeleton, and we then snap back the positions to the vertical planes afterwards.
  // Points of non-vertical faces that are shared with vertical faces are also moved.
  //
  // point id -> the vertical contour edges whose (lateral) faces contain this point
  CGAL::unordered_flat_map<std::size_t, std::vector<HDS_Halfedge_const_handle> > m_vertical_slabs;
#endif

  Slope m_slope;
  FT m_vertical_weight;
  std::optional<FT> m_height;

  std::vector<Point_3> m_points;
  std::vector<std::vector<std::size_t> > m_faces;
#ifdef CGAL_SLS_OUTPUT_FILES
  std::size_t m_current_cid;
  std::vector<std::size_t> m_colors;
#endif

public:
  Extrusion_builder(const SK& sk, const XK& xk, bool verbose = false)
    : m_sk(sk), m_xk(xk), m_verbose(verbose)
  { }

private:
  void reset_local_data()
  {
    m_slope = Slope::UNKNOWN;
    m_height = std::nullopt;
    m_vertex_point_ids.clear();
    m_halfedge_point_ids.clear();
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    m_vertical_slabs.clear();
#endif
    m_points.clear();
    m_faces.clear();
#ifdef CGAL_SLS_OUTPUT_FILES
    m_current_cid = 0;
    m_colors.clear();
#endif
  }

#ifdef CGAL_SLS_OUTPUT_FILES
  // Writes the soup with one color per triangulated border (bottom, top, lateral faces).
  // Must be called before `m_points` and `m_faces` are moved out.
  void write_colored_faces(const std::string& filename = "colored_faces.off") const
  {
    CGAL_precondition(m_colors.size() == m_faces.size());

    auto generate_colors = [](std::size_t n, float saturation = 0.70f, float lightness = 0.55f)
    {
      std::vector<std::array<int, 3>> palette;
      palette.reserve(n);

      // Start at a random hue so every run gives a fresh palette
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<float> dist(0.0f, 360.0f);
      float startHue = dist(gen);

      // Golden angle in degrees (~137.508) maximizes perceptual distance
      constexpr float goldenAngle = 137.507764f;

      for (std::size_t i = 0; i < n; ++i) {
        float h = std::fmod(startHue + float(i) * goldenAngle, 360.0f);

        // HSL to RGB conversion
        float c = (1.0f - std::abs(2.0f * lightness - 1.0f)) * saturation;
        float x = c * (1.0f - std::abs(std::fmod(h / 60.0f, 2.0f) - 1.0f));
        float m = lightness - c / 2.0f;

        float r = 0, g = 0, b = 0;
        if (h < 60)       { r = c; g = x; b = 0; }
        else if (h < 120) { r = x; g = c; b = 0; }
        else if (h < 180) { r = 0; g = c; b = x; }
        else if (h < 240) { r = 0; g = x; b = c; }
        else if (h < 300) { r = x; g = 0; b = c; }
        else              { r = c; g = 0; b = x; }

        palette.push_back({ int(std::round((r + m) * 255.0f)),
                            int(std::round((g + m) * 255.0f)),
                            int(std::round((b + m) * 255.0f)) });
      }
      return palette;
    };

    const auto palette = generate_colors(m_current_cid);

    std::ofstream outc(filename);
    outc.precision(17);
    outc << "COFF\n";
    outc << m_points.size() << " " << m_faces.size() << " 0\n";

    for (const Point_3& p : m_points)
      outc << p << " 255 255 255\n";

    for (std::size_t pos=0; pos<m_faces.size(); ++pos)
    {
      const std::vector<std::size_t>& f = m_faces[pos];
      const auto& c = palette[m_colors[pos]];
      outc << f.size();
      for (std::size_t vid : f)
        outc << " " << vid;
      outc << " " << c[0] << " " << c[1] << " " << c[2] << "\n";
    }
  }
#endif

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

public:
  template <typename K> // XK or EK
  typename K::Point_3 compute_skeleton_vertex_point(HDS_Vertex_const_handle vertex) const
  {
    const double s = (m_height.has_value() && is_negative(*m_height)) ? -1. : 1.;

    if constexpr(std::is_same_v<K, XK>)
    {
      FT time = vertex->time(); // FT here is XK::FT
      if (m_height.has_value())
        time = boost::algorithm::clamp(time, FT(0), CGAL::abs(*m_height));
      return { vertex->point().x(), vertex->point().y(), s * time };
    }
    else
    {
      // we requested the exact skeleton vertex position, so we need to convert the trisegment
      // and recompute the exact position.
      static_assert(is_exact_kernel<K>);
      using K_straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
      using Items_converter = CGAL::Straight_skeleton_items_converter_2<Straight_skeleton_2, K_straight_skeleton_2>;
      const auto k_vertex = Items_converter()(vertex);
      const typename K::FT& time = k_vertex.time();
      // if K is not XK, K is exact, so the time is correct
      CGAL_assertion(!m_height.has_value() || time <= CGAL::abs(*m_height));
      return { k_vertex.point().x(), k_vertex.point().y(), s * time };
    }
  }

  // @todo there is no caching of EK points at all at the moment...
  template <typename K> // XK or EK
  std::pair<typename K::Point_3, std::size_t>
  skeleton_vertex_point_and_id(HDS_Vertex_const_handle vertex)
  {
    CGAL::Cartesian_converter<K, XK> to_xk;

    const std::size_t n = m_points.size();
    auto res = m_vertex_point_ids.emplace(vertex, n);
    if(res.second) // successful insertion, never investigated this vertex before
    {
      const auto& kp = compute_skeleton_vertex_point<K>(vertex);

      // `m_points` is always XK.
      // In fact, `m_points` is always SK, but in the conversion EK->XK->SK, one of the cast is trivial
      // because either we use EK everywhere and XK == EK, or we only use it locally, and SK == XK
      m_points.push_back(to_xk(kp));

      return { kp, n };
    }

    // here and below, the insertion was *not* successful, meaning m_points already has that vertex
    const std::size_t id = res.first->second;

    if constexpr(is_exact_kernel<K>)
    {
      // if we have gone through the trouble of using an exact kernel, the position is likely better
      // than the one computed with XK, so update the XK position with the (converted-back) value.
      const auto& kp = compute_skeleton_vertex_point<K>(vertex);
      m_points[id] = to_xk(kp);

      return { kp, id };
    }
    else
    {
      return { m_points[id], id };
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////

  template <typename K> // XK or EK
  typename K::Point_3 compute_skeleton_halfedge_point(HDS_Halfedge_const_handle halfedge,
                                                      const FT& offset) const
  {
    using K_Point_2 = typename K::Point_2;
    using K_OB_traits = CGAL::Polygon_offset_builder_traits_2<K>;
    using K_Segment_2 = typename K_OB_traits::Segment_2_with_ID;
    using K_straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
    using Items_converter = CGAL::Straight_skeleton_items_converter_2<Straight_skeleton_2, K_straight_skeleton_2>;

    CGAL::Cartesian_converter<XK, K> to_k;
    auto make_segment = [&](HDS_Halfedge_const_handle contour_edge)
    {
      const K_Point_2& source = to_k(contour_edge->opposite()->vertex()->point());
      const K_Point_2& target = to_k(contour_edge->vertex()->point());
      const auto segment = K().construct_segment_2_object()(source, target);
      return K_Segment_2(segment, contour_edge->id());
    };

    HDS_Halfedge_const_handle contour_a = halfedge->defining_contour_edge();
    HDS_Halfedge_const_handle contour_b = halfedge->opposite()->defining_contour_edge();

    typename K_OB_traits::Trisegment_2_ptr seed_event;
    if(halfedge->is_inner_bisector())
    {
      HDS_Vertex_const_handle seed = (halfedge->slope() == CGAL::POSITIVE) ? halfedge->opposite()->vertex()
                                                                           : halfedge->vertex();

      // The trisegment is only read here, so there is no need to (deep) copy it if K == XK
      if constexpr(std::is_same_v<K, XK>)
        seed_event = seed->trisegment();
      else
        seed_event = Items_converter()(seed->trisegment());
    }

    K_OB_traits ob_traits;
    auto point = ob_traits.construct_offset_point_2_object()(to_k(CGAL::abs(offset)),
                                                             make_segment(contour_a),
                                                             to_k(contour_a->weight()),
                                                             make_segment(contour_b),
                                                             to_k(contour_b->weight()),
                                                             seed_event);
    CGAL_postcondition(point.has_value());

    return { point->x(), point->y(), offset };
  }

  template <typename K> // XK or EK
  std::pair<typename K::Point_3, std::size_t>
  skeleton_halfedge_point_and_id(HDS_Halfedge_const_handle halfedge)
  {
    CGAL::Cartesian_converter<K, XK> to_xk;

    HDS_Halfedge_const_handle canonical = (halfedge < halfedge->opposite()) ? halfedge
                                                                            : halfedge->opposite();

    const std::size_t n = m_points.size();
    auto res = m_halfedge_point_ids.emplace(canonical, n);
    if(res.second) // successful insertion, never investigated this halfedge before
    {
      const auto& kp = compute_skeleton_halfedge_point<K>(canonical, *m_height);
      m_points.push_back(to_xk(kp));
      return { kp, n };
    }

    // here and below, the insertion was *not* successful, meaning m_points already has that vertex
    const std::size_t id = res.first->second;

    if constexpr(is_exact_kernel<K>)
    {
      // if we have gone through the trouble of using an exact kernel, the position is likely better
      // than the one computed with XK, so update the XK position with the (converted-back) value.
      const auto& kp = compute_skeleton_halfedge_point<K>(canonical, *m_height);
      m_points[id] = to_xk(kp);

      return { kp, id };
    }
    else
    {
      return { m_points[id], id };
    }
  }

public:
  // A border is a set of directed edges, given as pairs of (point, point id).
  // The edges must form closed cycles, but the cycles do not need to be given in any particular order:
  // this allows to describe a domain with holes, or with parts that touch at a vertex.
  template <typename K>
  using Border_point = std::pair<typename K::Point_3, std::size_t>;

  template <typename K>
  using Border_edge = std::pair<Border_point<K>, Border_point<K> >;

  template <typename K>
  using Border_edges = std::vector<Border_edge<K> >;

  template <typename K>
  static void add_border_loop(const std::vector<Border_point<K> >& loop,
                              Border_edges<K>& edges)
  {
    for(std::size_t i=0; i<loop.size(); ++i)
      edges.emplace_back(loop[i], loop[(i+1) % loop.size()]);
  }

private:
  // Check that every border edge is an edge of a triangle, in the expected direction (this also rejects
  // a border edge that the CDT would have split at a vertex of another border).
  // The faces of the CDT are counterclockwise when seen from above: they go along border edges
  // in the same direction as the borders if the domain is on the left of the borders.
  template <typename K>
  static bool are_border_edges_oriented_triangle_edges(const std::vector<std::array<std::size_t, 3> >& triangles,
                                                        const Border_edges<K>& border_edges,
                                                        const bool same_direction)
  {
    // edges as (oriented) pairs of point ids
    using Id_edge = std::pair<std::size_t, std::size_t>;

    CGAL::unordered_flat_set<Id_edge, boost::hash<Id_edge> > triangle_halfedges;
    for(const auto& t : triangles)
      for(int i=0; i<3; ++i)
        triangle_halfedges.emplace(t[i], t[(i+1)%3]);

    for(const Border_edge<K>& e : border_edges)
    {
      const std::size_t source_id = e.first.second;
      const std::size_t target_id = e.second.second;
      const Id_edge h = same_direction ? Id_edge(source_id, target_id) : Id_edge(target_id, source_id);
      if(triangle_halfedges.count(h) == 0)
        return false;
    }

    return true;
  }

public:
  // Triangulates the domain bounded by the borders, and adds the triangles to the output faces.
  // Returns `false` if the triangulation fails, which  only happen with an inexact kernel.
  template <typename K> // XK or EK
  bool triangulate_borders(// m_points contain the XK embedding, but if XK != EK, we want the EK embedding
                           const Border_edges<K>& border_edges,
                           const bool invert_faces,
                           // whether the domain is on the left of the borders when seen from above
                           const bool domain_on_left)
  {
    // Projection_traits_xy_3 has benefits: it's an exact projection so we don't worse potentially
    // already inexact constructions and we avoid issues of consistency between faces for artificial
    // vertices. There are also downsides: the triangles in the output are not necessarily nicely
    // shaped, yielding higher chances of embedding issues at the very end (XK -> IK).
    using Itag = CGAL::No_constraint_intersection_tag;
    using PK = CGAL::Projection_traits_xy_3<K>;
    using Vbb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, PK>;
    using Vb = CGAL::Triangulation_vertex_base_2<PK, Vbb>;
    using Fb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using CDT = CGAL::Constrained_Delaunay_triangulation_2<PK, TDS, Itag>;
    using Vertex_handle = typename CDT::Vertex_handle;
    using Face_handle = typename CDT::Face_handle;

    CDT cdt;

    CGAL::unordered_flat_map<std::size_t, Vertex_handle> id_vertices;
    auto insert = [&](const Border_point<K>& bp) -> Vertex_handle
    {
      auto res = id_vertices.emplace(bp.second, Vertex_handle());
      if(res.second) // first time that this point is seen
      {
        res.first->second = cdt.insert(bp.first);
        res.first->second->info() = bp.second;
      }
      return res.first->second;
    };

    try
    {
      for(const Border_edge<K>& e : border_edges)
      {
        // two statements, so that the insertion order does not depend on the compiler
        const Vertex_handle va = insert(e.first);
        const Vertex_handle vb = insert(e.second);
        cdt.insert_constraint(va, vb);
      }
    }
    catch(const typename CDT::Intersection_of_constraints_exception&)
    {
      if(m_verbose)
        std::cerr << "Warning: the borders intersect in the triangulation (K = " << typeid(K).name() << ")"
                  << std::endl;
      return false;
    }

    // With an inexact kernel, distinct points can have the same projection, and the CDT then merges them
    // into a single vertex (pinching the domain, or collapsing an edge, whose constraint is then ignored):
    // every point must have its own vertex
    if(cdt.number_of_vertices() != id_vertices.size())
    {
      if(m_verbose)
        std::cerr << "Warning: distinct points of the borders are merged in the triangulation" << std::endl;
      return false;
    }

    // CDT has been built without issues, extract the faces in domain
#ifdef CGAL_SLS_DEBUG_DRAW
    // CGAL::draw(cdt);
#endif

    CGAL::unordered_flat_map<Face_handle, bool> in_domain_map;
    boost::associative_property_map<CGAL::unordered_flat_map<Face_handle, bool> > in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(cdt, in_domain);

    std::vector<std::array<std::size_t, 3> > triangles;
    for(Face_handle f : cdt.finite_face_handles())
    {
      if(!get(in_domain, f))
        continue;

      if(invert_faces)
        triangles.push_back({f->vertex(0)->info(), f->vertex(2)->info(), f->vertex(1)->info()});
      else
        triangles.push_back({f->vertex(0)->info(), f->vertex(1)->info(), f->vertex(2)->info()});
    }

    // With an inexact kernel, a sliver of the domain can be inverted (e.g., a top face of almost zero width,
    // when the offset is almost at the time of a horizontal ridge): the CDT then gives a valid triangulation
    // of the borders, but whose orientation is not consistent with the faces that share the borders.
    // Check that every border edge is an edge of a triangle, in the expected direction (this also rejects
    // a border edge that the CDT would have split at a vertex of another border).
    // The faces of the CDT are counterclockwise when seen from above: they go along border edges
    // in the same direction as the borders if the domain is on the left of the borders.
    if(!are_border_edges_oriented_triangle_edges<K>(triangles, border_edges, domain_on_left != invert_faces))
    {
      if(m_verbose)
        std::cerr << "Warning: the triangulation is not consistent with the orientation of the borders" << std::endl;
      return false;
    }

    for(const auto& t : triangles)
    {
      m_faces.push_back({t[0], t[1], t[2]});
#ifdef CGAL_SLS_OUTPUT_FILES
      m_colors.push_back(m_current_cid);
#endif
    }

#ifdef CGAL_SLS_OUTPUT_FILES
    ++m_current_cid;
#endif

    return true;
  }

  void collect_contour_halfedges(const Straight_skeleton_2& ss,
                                 std::vector<std::vector<HDS_Halfedge_const_handle> >& chs)
  {
    CGAL::unordered_flat_map<HDS_Halfedge_const_handle, bool> visited;
    for(const HDS_Halfedge_const_handle hds_h : CGAL::halfedges(ss))
    {
      if(visited.find(hds_h) != visited.end() || !hds_h->is_border())
        continue;

      std::vector<HDS_Halfedge_const_handle> contour;
      HDS_Halfedge_const_handle hds_h0 = hds_h;
      do
      {
        contour.push_back(hds_h0);
        visited[hds_h0] = true;
        hds_h0 = hds_h0->next();
      }
      while(hds_h0 != hds_h);

      chs.push_back(std::move(contour));
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////

  bool construct_bot_horizontal_faces(const std::vector<std::vector<HDS_Halfedge_const_handle> >& chs)
  {
    // convert from the contours to a vector of points and a vector of point ids
    Border_edges<XK> border_edges;
    for(const std::vector<HDS_Halfedge_const_handle>& ch : chs)
    {
      std::vector<Border_point<XK> > loop;
      for(HDS_Halfedge_const_handle hds_h : ch)
        loop.push_back(skeleton_vertex_point_and_id<XK>(hds_h->vertex()));
      add_border_loop<XK>(loop, border_edges);
    }

    // The faces are built in the plane z=0 with upwards normal.
    // The borders are the border halfedges of the skeleton(s), on the other side of the skeleton faces:
    // the bottom face is on their right for an interior skeleton (inward), and on their left for
    // exterior skeletons (outward).
    const bool invert_faces = !m_height.has_value() || is_positive(*m_height);
    const bool domain_on_left = (m_slope == Slope::OUTWARD);
    return triangulate_borders<XK>(border_edges, invert_faces, domain_on_left);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////

  template <typename K> // XK or EK
  bool triangulate_top_contours_in_kernel(const std::vector<std::vector<HDS_Halfedge_const_handle> >& chs)
  {
    // convert from the contours to border edges
    Border_edges<K> border_edges;
    for(const std::vector<HDS_Halfedge_const_handle>& ch : chs)
    {
      std::vector<Border_point<K> > loop;
      for(HDS_Halfedge_const_handle hds_h : ch)
        loop.push_back(skeleton_halfedge_point_and_id<K>(hds_h));
      add_border_loop<K>(loop, border_edges);
    }

    // The top faces are on the left of the offset contours of an interior skeleton (inward),
    // and on their right for exterior skeletons (outward).
    const bool invert_faces = m_height.has_value() && is_negative(*m_height);
    const bool domain_on_left = (m_slope != Slope::OUTWARD);
    return triangulate_borders<K>(border_edges, invert_faces, domain_on_left);
  }

  // this dispatches to the appropriate kernel
  bool triangulate_top_contours(const std::vector<std::vector<HDS_Halfedge_const_handle> >& chs,
                                const bool exact)
  {
    if constexpr(std::is_floating_point<FT>::value)
    {
      if(exact)
      {
        using EK = typename Extrusion_kernel_selector<XK>::type;
        return triangulate_top_contours_in_kernel<EK>(chs);
      }
    }

    return triangulate_top_contours_in_kernel<XK>(chs);
  }

  // this is the robust call: try with XK, and if it fails, use EK
  bool construct_top_horizontal_faces(const std::vector<std::vector<HDS_Halfedge_const_handle> >& chs)
  {
    if (chs.empty())
      return true;

    bool success;
    success = triangulate_top_contours(chs, false /*exact*/);
    if constexpr(std::is_floating_point<FT>::value)
    {
      if(!success)
        success = triangulate_top_contours(chs, true /*exact*/);
    }

    CGAL_postcondition(success);
    if(!success && m_verbose)
      std::cerr << "Error: Failed to triangulate skeleton top" << std::endl;

    return success;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////

  static bool is_artificial_bisector(HDS_Halfedge_const_handle h)
  {
    return h->face() == h->opposite()->face();
  }

  // @todo The crossings of the offset with artificial bisectors are in the interior of offset segments,
  // and could be ignored, both here and in top faces (see `Skeleton_offset_correspondence_builder_visitor`).
  // With an inexact kernel, this makes the output a bit less likely to self-intersect when the offset is
  // (almost) at the time of an endpoint of an artificial bisector.
  template <typename K> // XK or EK
  bool triangulate_skeleton_face_in_kernel(HDS_Face_handle hds_f,
                                           const Offset_builder_no_visitor& ob,
                                           const bool invert_faces)
  {
    HDS_Halfedge_const_handle hds_h = hds_f->halfedge(), done = hds_h;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    HDS_Halfedge_const_handle contour_h = hds_h->defining_contour_edge();
    CGAL_assertion(hds_h == contour_h);
    const bool is_vertical = (contour_h->weight() == m_vertical_weight);
#endif

    // Note that this comparison uses XK (and not K), but XK is supposed to have exact predicates.
    //
    // @todo we probably don't even need to do this check here again:
    // if we have contoured at z=height, we already know the relevant halfedges
    auto compare_time_to_offset = [&](HDS_Vertex_const_handle node) -> CGAL::Comparison_result
    {
      if(!m_height.has_value())
        return CGAL::EQUAL; // meaningless, but just to avoid getting a point on the halfedge
      else if(node->is_contour())
        return CGAL::LARGER; // abs_height is > 0 and contour is 0
      else
        return ob.Compare_offset_against_event_time(CGAL::abs(*m_height), node);
    };

    // The points of the boundary cycle, each with whether the edge from the previous point
    // to this point is an edge of the border (that is, not a part of an artificial bisector)
    std::vector<std::pair<Border_point<K>, bool> > points;

    do
    {
      HDS_Vertex_const_handle hds_sv = hds_h->opposite()->vertex();
      HDS_Vertex_const_handle hds_tv = hds_h->vertex();

      // abs(height) compared to node->time
      // smaller means the node is farther than the offset
      const CGAL::Comparison_result sc = compare_time_to_offset(hds_sv);
      const CGAL::Comparison_result tc = compare_time_to_offset(hds_tv);

      // we skip artificial bisectors because CDT2 doesn't like duplicate constraints
      const bool is_artificial = is_artificial_bisector(hds_h);

      if(sc != tc && sc != CGAL::EQUAL && tc != CGAL::EQUAL) // intersect the offset in the interior
      {
        // hds_h gets canonicalized in the point/id getter
        const Border_point<K> bp = skeleton_halfedge_point_and_id<K>(hds_h);

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
        if(is_vertical)
          add_vertical_slab(bp.second, contour_h);
#endif

        // Going down (the source is above the offset): the edge from the previous point, where the cycle
        // went above the offset, is an offset segment. Going up: it is the part of `hds_h` from its source.
        const bool is_offset_segment = (sc == CGAL::SMALLER);
        points.emplace_back(bp, is_offset_segment || !is_artificial);
      }

      if(tc != CGAL::SMALLER)
      {
        const Border_point<K> bp = skeleton_vertex_point_and_id<K>(hds_tv);

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
        // contour vertices are already on the vertical planes
        if(is_vertical && !hds_tv->is_contour())
          add_vertical_slab(bp.second, contour_h);
#endif

        // The edge from the previous point is an offset segment if `hds_h` comes down from above the offset
        // directly to a target that is on the offset (no crossing was emitted before it). Otherwise, it is
        // the part of `hds_h` from its source or from the crossing.
        const bool is_offset_segment = (sc == CGAL::SMALLER && tc == CGAL::EQUAL);
        points.emplace_back(bp, is_offset_segment || !is_artificial);
      }

      hds_h = hds_h->next();
    }
    while(hds_h != done);

    Border_edges<K> border_edges;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Border_point<K>& prev_point = points[i].first;
      const auto& [point, is_border_edge] = points[(i + 1) % points.size()];

      // e.g., an offset segment from a vertex that is on the offset to itself, when the boundary cycle
      // goes up from the vertex to a hole that is above the offset except at that vertex, and back
      if(is_border_edge && prev_point.second != point.second)
        border_edges.emplace_back(prev_point, point);
    }

    if(border_edges.size() < 3)
    {
      if(m_verbose)
        std::cerr << "Warning: sm_vs has size 1 or 2: offset crossing face at a single point?" << std::endl;
      return true;
    }

    // the borders follow the boundary cycle of the skeleton face, which is on their left
    return triangulate_borders<K>(border_edges, invert_faces, true /*domain on the left*/);
  }

  // this dispatches to the appropriate kernel
  bool triangulate_skeleton_face(HDS_Face_handle hds_f,
                                 const Offset_builder_no_visitor& ob,
                                 const bool invert_faces,
                                 const bool use_exact)
  {
    if constexpr(!is_exact_kernel<XK>)
    {
      if(use_exact)
      {
        using EK = typename Extrusion_kernel_selector<XK>::type;
        return triangulate_skeleton_face_in_kernel<EK>(hds_f, ob, invert_faces);
      }
    }

    return triangulate_skeleton_face_in_kernel<XK>(hds_f, ob, invert_faces);
  }

  // this is the robust call: try with XK, and if it fails, use EK
  void construct_lateral_faces(const Straight_skeleton_2& ss,
                               const bool ignore_frame_faces = false,
                               const bool invert_faces = false)
  {
    Offset_builder_no_visitor ob(ss);

    const HDS& hds = static_cast<const HDS&>(ss);

    std::size_t fc = 0;
    for(const HDS_Face_handle hds_f : CGAL::faces(hds))
    {
      // If they exist (exterior skeleton), the first four faces of the SLS correspond
      // to the outer frame, and should be ignored.
      if(ignore_frame_faces && fc++ < 4)
        continue;

      bool success;
      success = triangulate_skeleton_face(hds_f, ob, invert_faces, false /*exact*/);
      if constexpr(std::is_floating_point<FT>::value)
      {
        if(!success)
          success = triangulate_skeleton_face(hds_f, ob, invert_faces, true /*exact*/);
      }

      CGAL_postcondition(success);
      if(!success && m_verbose)
        std::cerr << "Error: Failed to triangulate lateral face" << std::endl;
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
  void add_vertical_slab(const std::size_t pid,
                         HDS_Halfedge_const_handle contour_h)
  {
    std::vector<HDS_Halfedge_const_handle>& slabs = m_vertical_slabs[pid];
    if(std::find(slabs.begin(), slabs.end(), contour_h) == slabs.end())
      slabs.push_back(contour_h);
  }

  // orthogonal projection onto the supporting line of the (vertical) contour edge `slab`
  Point_2 project_onto_slab(const Point_2& p,
                            HDS_Halfedge_const_handle slab) const
  {
    const Point_2& s = slab->opposite()->vertex()->point();
    const Point_2& t = slab->vertex()->point();

    // axis-aligned slabs: keep the exact coordinate of the input
    if(s.x() == t.x())
      return { s.x(), p.y() };
    if(s.y() == t.y())
      return { p.x(), s.y() };

    return m_xk.construct_projected_point_2_object()(m_xk.construct_line_2_object()(s, t), p);
  }

  // Moves the points of the vertical faces onto the vertical planes above their contour edges.
  void snap_to_vertical_slabs()
  {
    auto cl2 = m_xk.construct_line_2_object();
    auto collinear = m_xk.collinear_2_object();
    auto cp2 = m_xk.construct_point_2_object();
    auto cp3 = m_xk.construct_point_3_object();
    auto inter_2 = m_xk.intersect_2_object();

    auto have_same_supporting_line = [&](HDS_Halfedge_const_handle h0,
                                         HDS_Halfedge_const_handle h1)
    {
      const Point_2& s0 = h0->opposite()->vertex()->point();
      const Point_2& t0 = h0->vertex()->point();
      return collinear(s0, t0, h1->opposite()->vertex()->point()) &&
             collinear(s0, t0, h1->vertex()->point());
    };

    for(const auto& [pid, slabs] : m_vertical_slabs)
    {
      Point_3& p = m_points[pid];
      const Point_2 p2 = cp2(p.x(), p.y());

      if(slabs.size() > 2)
      {
        if(m_verbose)
          std::cerr << "Warning: point " << p << " is on " << slabs.size() << " distinct vertical slabs??" << std::endl;
      }

      // keep a single slab per supporting line (e.g., consecutive collinear vertical edges)
      std::vector<HDS_Halfedge_const_handle> lines;
      for(HDS_Halfedge_const_handle w : slabs)
      {
        bool is_new = std::none_of(lines.begin(), lines.end(),
                                   [&](HDS_Halfedge_const_handle l) { return have_same_supporting_line(l, w); });
        if(is_new)
          lines.push_back(w);
      }

      std::optional<Point_2> snapped;
      if(lines.size() == 1)
      {
        snapped = project_onto_slab(p2, lines[0]);
      }
      else if(lines.size() >= 2)
      {
        // On paper, one could build a polygon very large such that we get a node
        // with 3 edges assigned the fake vertical weight...
        // @fixme scale the weight by the characteristic size of the polgyon?
        CGAL_warning_msg(lines.size() == 2, "a point is on more than two vertical slabs");

        // check the simple case of consecutive edges to avoid computations
        if(lines[0]->vertex() == lines[1]->opposite()->vertex())
        {
          snapped = lines[0]->vertex()->point();
        }
        else if(lines[1]->vertex() == lines[0]->opposite()->vertex())
        {
          snapped = lines[1]->vertex()->point();
        }
        else
        {
          const Line_2 l0 = cl2(lines[0]->opposite()->vertex()->point(), lines[0]->vertex()->point());
          const Line_2 l1 = cl2(lines[1]->opposite()->vertex()->point(), lines[1]->vertex()->point());
          const auto inter = inter_2(l0, l1);
          // The lines are distinct, but they can be parallel (same reason as above)
          CGAL_warning_msg(bool(inter), "parallel (but distinct) lines?...");
          if(inter)
          {
            if(const Point_2* ip = std::get_if<Point_2>(&*inter))
              snapped = *ip;
          }
        }
      }

      if(snapped)
        p = cp3(snapped->x(), snapped->y(), p.z());
    }
  }
#endif // CGAL_SLS_SNAP_TO_VERTICAL_SLABS

public:
  // This is roughly "CGAL::create_interior_weighted_skeleton_and_offset_polygons_with_holes_2()",
  // but we want to know the intermediate straight skeleton to build the lateral faces of the 3D mesh.
  template <typename PolygonWithHoles, typename WeightRange, typename IFT,
            typename PointRange, typename FaceRange>
  bool inward_construction(const PolygonWithHoles& pwh,
                           const WeightRange& speeds,
                           const IFT i_vertical_weight,
                           const SFT s_height,
                           PointRange& points,
                           FaceRange& faces)
  {
    using IPoint_2 = typename boost::range_value<typename PolygonWithHoles::Polygon_2>::type;
    using IK = typename Kernel_traits<IPoint_2>::Kernel;

#ifdef CGAL_SLS_DEBUG_DRAW
    CGAL::draw(pwh);
#endif

    reset_local_data();

    CGAL::Cartesian_converter<IK, SK> to_sk;
    CGAL::Cartesian_converter<SK, XK> to_xk;

    m_slope = Slope::INWARD;
    m_vertical_weight = to_xk(to_sk(i_vertical_weight));
    if(!is_zero(s_height))
      m_height = to_xk(s_height);

    const bool is_default_extrusion_height = !(m_height.has_value());
    const bool is_partial_skeleton = (!is_default_extrusion_height && pwh.number_of_holes() == 0);

    // @partial_wsls_pwh interior SLS of weighted polygons with holes can have skeleton faces with holes
    // The current postprocessing is in the function EnforceSimpleConnectedness, but it has not yet
    // been made compatible with partial skeletons.
    S_Straight_skeleton_2_ptr s_ss_ptr;

#ifdef CGAL_SLS_EXTRUSION_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    if(!is_partial_skeleton)
    {
      s_ss_ptr = CGAL::create_interior_weighted_straight_skeleton_2(
                   CGAL_SS_i::vertices_begin(pwh.outer_boundary()),
                   CGAL_SS_i::vertices_end(pwh.outer_boundary()),
                   pwh.holes_begin(), pwh.holes_end(),
                   std::begin(speeds[0]), std::end(speeds[0]),
                   std::next(std::begin(speeds)), std::end(speeds),
                   m_sk);
    }
    else
    {
      s_ss_ptr = CGAL_SS_i::create_partial_interior_weighted_straight_skeleton_2(
                   CGAL::abs(s_height),
                   CGAL_SS_i::vertices_begin(pwh.outer_boundary()),
                   CGAL_SS_i::vertices_end(pwh.outer_boundary()),
                   pwh.holes_begin(), pwh.holes_end(),
                   std::begin(speeds[0]), std::end(speeds[0]),
                   std::next(std::begin(speeds)), std::end(speeds),
                   m_sk);
    }

#ifdef CGAL_SLS_EXTRUSION_TIMERS
    timer.stop();
    std::cout << "Straight skeleton construction took " << timer.time() << " seconds" << std::endl;
#endif

    if(!s_ss_ptr)
    {
      if(m_verbose)
        std::cerr << "Error: encountered an error during skeleton construction" << std::endl;
      return false;
    }

#ifdef CGAL_SLS_DEBUG_DRAW
    // Straight_skeletons_2::IO::print_straight_skeleton(*s_ss_ptr);
    CGAL::draw(*s_ss_ptr);
#endif

#ifdef CGAL_SLS_EXTRUSION_TIMERS
    timer.reset();
    timer.start();
#endif

    CGAL_assertion(s_ss_ptr->is_valid(is_partial_skeleton));
    Straight_skeleton_2_ptr ss_ptr = CGAL::convert_straight_skeleton_2<Straight_skeleton_2>(s_ss_ptr);
    CGAL_assertion(ss_ptr && ss_ptr->is_valid(is_partial_skeleton));

#ifdef CGAL_SLS_EXTRUSION_TIMERS
    timer.stop();
    std::cout << "Conversion took " << timer.time() << " seconds" << std::endl;
    timer.reset();
    timer.start();
#endif

    // We construct face by face independently because that enables us to switch to an exact kernel
    // only on-demand as opposed to recomputing the whole skeleton with an exact construction kernel
    // as soon as we get a single, local issue.

    // bottom faces (z=0)
    {
      // grab the outer boundary contour and the hole contours
      std::vector<std::vector<HDS_Halfedge_const_handle> > contour_halfedges;
      collect_contour_halfedges(*ss_ptr, contour_halfedges);
      CGAL_assertion(contour_halfedges.size() == 1 + pwh.number_of_holes());
      construct_bot_horizontal_faces(contour_halfedges);
    }

    // top faces
    if(!is_default_extrusion_height)
    {
      std::vector<std::vector<HDS_Halfedge_const_handle> > contour_halfedges;
      Visitor visitor(*m_height, m_points, m_vertex_point_ids, m_halfedge_point_ids, contour_halfedges);
      Offset_builder ob(*ss_ptr, Offset_builder_traits(), visitor);
      ob.construct_offset_contours(CGAL::abs(*m_height), CGAL::Emptyset_iterator());
      construct_top_horizontal_faces(contour_halfedges);
    }

    // lateral faces
    {
      const bool invert_faces = m_height.has_value() && is_negative(*m_height);
      construct_lateral_faces(*ss_ptr, false /*no frame faces*/, invert_faces);
    }

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    snap_to_vertical_slabs();
#endif

#ifdef CGAL_SLS_EXTRUSION_TIMERS
    timer.stop();
    std::cout << "Offset construction took " << timer.time() << " seconds" << std::endl;
#endif

#ifdef CGAL_SLS_OUTPUT_FILES
    write_colored_faces();
#endif

    points = std::move(m_points);
    faces = std::move(m_faces);

    return true;
  }

  template <typename PolygonWithHoles, typename WeightRange, typename IFT,
            typename PointRange, typename FaceRange>
  bool outward_construction(const PolygonWithHoles& pwh,
                            const WeightRange& speeds,
                            const IFT i_vertical_weight,
                            const SFT s_height,
                            PointRange& points,
                            FaceRange& faces)
  {
    CGAL_precondition(!is_zero(s_height));

#ifdef CGAL_SLS_DEBUG_DRAW
    CGAL::draw(pwh);
#endif

    reset_local_data();

    using IPoint_2 = typename boost::range_value<typename PolygonWithHoles::Polygon_2>::type;
    using IK = typename Kernel_traits<IPoint_2>::Kernel;

    CGAL::Cartesian_converter<IK, SK> to_sk;
    CGAL::Cartesian_converter<SK, XK> to_xk;

    m_slope = Slope::OUTWARD;
    m_vertical_weight = to_xk(to_sk(i_vertical_weight));
    m_height = to_xk(s_height);

    std::vector<std::vector<HDS_Halfedge_const_handle> > bot_contour_halfedges;
    std::vector<std::vector<HDS_Halfedge_const_handle> > top_contour_halfedges;

    // The exterior of a polygon with holes is the exterior of its outer boundary,
    // and the interior of its inverted holes.

    // Start with the outer boundary
    Straight_skeleton_2_ptr ss_ptr; // here because we want to keep the HDS alive
    {
      S_Straight_skeleton_2_ptr s_ss_ptr = CGAL_SS_i::create_partial_exterior_weighted_straight_skeleton_2(
                                             CGAL::abs(s_height),
                                             CGAL_SS_i::vertices_begin(pwh.outer_boundary()),
                                             CGAL_SS_i::vertices_end(pwh.outer_boundary()),
                                             std::begin(speeds[0]), std::end(speeds[0]),
                                             m_sk);

      if(!s_ss_ptr)
      {
        if(m_verbose)
          std::cerr << "Error: encountered an error during outer skeleton construction" << std::endl;
        return false;
      }

#ifdef CGAL_SLS_DEBUG_DRAW
      Straight_skeletons_2::IO::print_straight_skeleton(*s_ss_ptr);
      CGAL::draw(*s_ss_ptr);
#endif

      CGAL_assertion(s_ss_ptr->is_valid(true /*partial*/));
      ss_ptr = CGAL::convert_straight_skeleton_2<Straight_skeleton_2>(s_ss_ptr);
      CGAL_assertion(ss_ptr && ss_ptr->is_valid(true /*partial*/));

      // accumulate contours for the bottom faces
      collect_contour_halfedges(*ss_ptr, bot_contour_halfedges);
      CGAL_assertion(bot_contour_halfedges.size() == 2);

      // filter the outer frame
      std::swap(bot_contour_halfedges[0], bot_contour_halfedges.back());
      bot_contour_halfedges.pop_back();

      // accumulate the result without triangulating anything yet: we need all the holes in one place
      Visitor visitor(*m_height, m_points, m_vertex_point_ids, m_halfedge_point_ids, top_contour_halfedges, true /*ignore the first contour*/);
      Offset_builder ob(*ss_ptr, Offset_builder_traits(), visitor);
      ob.construct_offset_contours(CGAL::abs(*m_height), CGAL::Emptyset_iterator());

      // nothing to purge from the top contour since we ignore it in the callback

      // lateral faces
      const bool invert_faces = is_positive(*m_height);
      construct_lateral_faces(*ss_ptr, true /*ignore frame faces*/, invert_faces);
    }

    // now, deal with the holes
    std::vector<Straight_skeleton_2_ptr> hole_ss_ptrs; // keep the HDSs alive

    std::size_t hole_id = 0;
    for(auto hit=pwh.holes_begin(); hit!=pwh.holes_end(); ++hit, ++hole_id)
    {
      using I_Polygon_2 = CGAL::cpp20::remove_cvref_t<decltype(*pwh.holes_begin())>;

      // reverse the hole and its weights
      std::vector<IPoint_2> hole(std::cbegin(*hit), std::cend(*hit));
      std::reverse(hole.begin(), hole.end());

      auto hole_speeds = speeds[1 + hole_id];
      std::reverse(hole_speeds.begin(), hole_speeds.end());
      std::rotate(hole_speeds.begin(), std::prev(hole_speeds.end()), hole_speeds.end()); // one to the right

      // this is roughly "CGAL::create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2()",
      // but we need the intermediate straight skeleton to build the lateral faces

      std::vector<I_Polygon_2> no_holes;
      std::vector<std::vector<IFT> > no_hole_speeds;
      S_Straight_skeleton_2_ptr s_ss_ptr = CGAL_SS_i::create_partial_interior_weighted_straight_skeleton_2(
                                             CGAL::abs(s_height),
                                             CGAL_SS_i::vertices_begin(hole),
                                             CGAL_SS_i::vertices_end(hole),
                                             std::cbegin(no_holes),
                                             std::cend(no_holes),
                                             std::cbegin(hole_speeds),
                                             std::cend(hole_speeds),
                                             std::cbegin(no_hole_speeds),
                                             std::cend(no_hole_speeds),
                                             m_sk);

      if(!s_ss_ptr)
      {
        if(m_verbose)
          std::cerr << "Error: encountered an error during hole skeleton construction" << std::endl;
        return false;
      }

#ifdef CGAL_SLS_DEBUG_DRAW
      Straight_skeletons_2::IO::print_straight_skeleton(*s_ss_ptr);
      CGAL::draw(*s_ss_ptr);
#endif

      CGAL_assertion(s_ss_ptr->is_valid(true /*partial*/));
      Straight_skeleton_2_ptr ss_ptr = CGAL::convert_straight_skeleton_2<Straight_skeleton_2>(s_ss_ptr);
      CGAL_assertion(ss_ptr && ss_ptr->is_valid(true /*partial*/));
      hole_ss_ptrs.push_back(ss_ptr);

      collect_contour_halfedges(*ss_ptr, bot_contour_halfedges);

      // accumulate the result without triangulating anything yet: we need all the holes too
      Visitor visitor(*m_height, m_points, m_vertex_point_ids, m_halfedge_point_ids, top_contour_halfedges);
      Offset_builder ob(*ss_ptr, Offset_builder_traits(), visitor);
      ob.construct_offset_contours(CGAL::abs(*m_height), CGAL::Emptyset_iterator());

      // no frame contour to drop here: we used an interior skeleton for the reversed hole

      // lateral faces
      const bool invert_faces = is_positive(*m_height);
      construct_lateral_faces(*ss_ptr, false /*no outer frame*/, invert_faces);
    }

    // horizontal faces
    construct_bot_horizontal_faces(bot_contour_halfedges);
    construct_top_horizontal_faces(top_contour_halfedges);

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    snap_to_vertical_slabs();
#endif

#ifdef CGAL_SLS_OUTPUT_FILES
    write_colored_faces();
#endif

    points = std::move(m_points);
    faces = std::move(m_faces);

    return true;
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// convert angles (in degrees) to weights, and handle vertical angles
template <typename FT, typename AngleRange>
void convert_angles(AngleRange& angles)
{
  CGAL_precondition(!angles.empty());

  auto angle_to_weight = [](const FT& angle) -> FT
  {
    CGAL_precondition(0 < angle && angle < 180);

    // @todo should this be an epsilon around 90°? As theta goes to 90°, tan(theta) goes to infinity
    // and thus we could get numerical issues (overflows) if the kernel is not exact
    if(angle == 90)
      return 0;
    else
      return std::tan(CGAL::to_double(angle * CGAL_PI / 180));
  };

  for(auto& contour_angles : angles)
    for(FT& angle : contour_angles)
      angle = angle_to_weight(angle);
}

// handle vertical angles (inf speed)
// returns
// - whether the weights are positive or negative (inward / outward)
// - whether the input weights are valid
// - the weight of vertical slabs
template <typename FT, typename WeightRange>
std::tuple<Slope, bool, FT>
preprocess_weights(WeightRange& weights,
                   bool verbose = false)
{
  CGAL_precondition(!weights.empty());

  Slope slope = Slope::UNKNOWN;

  FT max_value = 0; // non-inf, maximum absolute value
  for(auto& contour_weights : weights)
  {
    for(FT& w : contour_weights)
    {
      // '0' means a vertical slab, aka 90° angle (see preprocess_angles())
      if(w == 0)
        continue;

      // determine whether weights indicate all inward or all outward
      if(slope == Slope::UNKNOWN)
      {
        // w is neither 0 nor inf here
        slope = (w > 0) ? Slope::INWARD : Slope::OUTWARD;
      }
      else if(slope == Slope::INWARD && w < 0)
      {
        if(verbose)
          std::cerr << "Error: mixing positive and negative weights is not yet supported" << std::endl;
        return {Slope::UNKNOWN, false, FT(-1)};
      }
      else if(slope == Slope::OUTWARD && w > 0)
      {
        if(verbose)
          std::cerr << "Error: mixing positive and negative weights is not yet supported" << std::endl;
        return {Slope::UNKNOWN, false, FT(-1)};
      }

      // if we are going outwards, it is just an interior skeleton with opposite weights
      w = CGAL::abs(w);
      if(w > max_value)
        max_value = w;
    }
  }

  if(slope == Slope::UNKNOWN)
  {
    if(verbose)
      std::cerr << "Warning: all edges vertical?" << std::endl;
    slope = Slope::VERTICAL;
  }

  // Take a weight which is a large % of the max value to ensure there's no ambiguity
  //
  // Since the max value might not be very close to 90°, take the max between of the large-% weight
  // and the weight corresponding to an angle of 89.9999999°
  const FT weight_of_89d9999999 = 572957787.3425436; // tan(89.9999999°)
  const FT scaled_max = (std::max)(weight_of_89d9999999, FT(1e3 * max_value));

  for(auto& contour_weights : weights)
  {
    for(FT& w : contour_weights)
    {
      if(w == FT(0))
        w = scaled_max;
    }
  }

  return {slope, true, scaled_max};
}

template <typename PolygonWithHoles,
          typename WeightRange,
          typename PolygonMesh,
          typename NamedParameters>
bool extrude_skeleton(const PolygonWithHoles& pwh,
                      WeightRange& weights,
                      PolygonMesh& out,
                      const NamedParameters& np)
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using I_Polygon_2 = typename PolygonWithHoles::Polygon_2;
  using I_Point = typename boost::range_value<I_Polygon_2>::type;
  using IK = typename Kernel_traits<I_Point>::type;
  using IFT = typename IK::FT;

  using SK = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                          NamedParameters,
                                                          IK>::type;
  using SFT = typename SK::FT;

  using XK = typename internal_np::Lookup_named_param_def<internal_np::extrusion_traits_t,
                                                          NamedParameters,
                                                          SK>::type;
  using Point_3 = typename XK::Point_3;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  const bool verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);
  SK sk = choose_parameter<SK>(get_parameter(np, CGAL::internal_np::geom_traits));
  XK xk = choose_parameter<XK>(get_parameter(np, internal_np::extrusion_traits));

  if(verbose)
    std::cout << "XK is " << typeid(XK).name() << std::endl;

  constexpr bool has_height = !is_default_parameter<NamedParameters, internal_np::maximum_height_t>::value;
  CGAL_precondition(!has_height || !is_zero(choose_parameter(get_parameter(np, internal_np::maximum_height), 0)));
  CGAL_USE(has_height);

  Slope slope;
  bool is_valid_input;
  IFT i_vertical_weight;
  std::tie(slope, is_valid_input, i_vertical_weight) = preprocess_weights<IFT>(weights, verbose);

  if(!is_valid_input)
  {
    if(verbose)
      std::cerr << "Error: invalid input weights" << std::endl;
    return false;
  }

  if(verbose)
  {
    switch(slope)
    {
      case Slope::UNKNOWN: std::cout << "Slope is UNKNOWN??" << std::endl; break;
      case Slope::INWARD: std::cout << "Slope is INWARD" << std::endl; break;
      case Slope::OUTWARD: std::cout << "Slope is OUTWARD" << std::endl; break;
      case Slope::VERTICAL: std::cout << "Slope is VERTICAL" << std::endl; break;
    }
  }

  // End of preprocessing, start the actual skeleton computation

  const SFT s_height = choose_parameter(get_parameter(np, internal_np::maximum_height), 0);
  if(slope != Slope::INWARD && s_height == 0)
  {
    if(verbose)
      std::cerr << "Error: a maximal height must be specified when using an outward (and vertical) slopes" << std::endl;
    return false;
  }

  // build a soup, to be converted to a mesh afterwards
  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > faces;
  points.reserve(2 * pwh.outer_boundary().size()); // just a reasonable guess
  faces.reserve(2 * pwh.outer_boundary().size() + 2*pwh.number_of_holes());

  Extrusion_builder<SK, XK> builder(sk, xk, verbose);
  bool res;
  if(slope != Slope::OUTWARD) // INWARD or VERTICAL
    res = builder.inward_construction(pwh, weights, i_vertical_weight, s_height, points, faces);
  else
    res = builder.outward_construction(pwh, weights, i_vertical_weight, s_height, points, faces);

  if(!res)
    return false;

#ifdef CGAL_SLS_OUTPUT_FILES
  // This soup provides one connected component per edge of the input polygon
  CGAL::IO::write_polygon_soup("extruded_skeleton_soup.off", points, faces, CGAL::parameters::stream_precision(17));
#endif

  using VPM = typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type;
  VPM out_vpm = choose_parameter(get_parameter(np, CGAL::internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, out));

  using O_Point_3 = typename boost::property_traits<VPM>::value_type;
  using OK = typename CGAL::Kernel_traits<O_Point_3>::type;

  CGAL::Cartesian_converter<XK, OK> to_ok;
  std::vector<O_Point_3> output_points;
  output_points.reserve(points.size());
  for(const Point_3& p : points)
    output_points.push_back(to_ok(p));

#ifdef CGAL_SLS_OUTPUT_FILES
  CGAL::IO::write_polygon_soup("result_soup.off", output_points, faces,
                               CGAL::parameters::stream_precision(17));
#endif

  CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(faces));
  PMP::polygon_soup_to_polygon_mesh(output_points, faces, out,
                                    CGAL::parameters::all_default(),
                                    CGAL::parameters::vertex_point_map(out_vpm));

#ifdef CGAL_SLS_DEBUG_DRAW
  CGAL::draw(out);
#endif

#ifdef CGAL_SLS_OUTPUT_FILES
  CGAL::IO::write_polygon_mesh("result.off", out, CGAL::parameters::stream_precision(17));
#endif

  CGAL_postcondition(is_valid_polygon_mesh(out));
  CGAL_postcondition(is_closed(out));

  return true;
}

} // namespace internal
} // namespace Straight_skeleton_extrusion

// Documented in the Straight_skeleton_2 package
template <typename PolygonWithHoles,
          typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
bool extrude_skeleton(const PolygonWithHoles& pwh,
                      PolygonMesh& out,
                      const NamedParameters& np = parameters::default_values(),
                      std::enable_if_t<CGAL_SS_i::has_Hole_const_iterator<PolygonWithHoles>::value>* = nullptr)
{
  namespace SSEI = Straight_skeleton_extrusion::internal;

  using Polygon_2 = typename PolygonWithHoles::General_polygon_2;
  using K = typename Kernel_traits<typename boost::range_value<Polygon_2>::type>::type;
  using FT = typename K::FT;
  using Default_speed_type = std::vector<std::vector<FT> >;

  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter_reference;

  const bool has_weights = !(is_default_parameter<NamedParameters, internal_np::weights_param_t>::value);
  const bool has_angles = !(is_default_parameter<NamedParameters, internal_np::angles_param_t>::value);

  Default_speed_type def_speeds; // never used, but needed for compilation

  if(has_weights)
  {
    // not using "Angles" here is on purpose, I want to copy the range but the NP is ref only
    auto weights = choose_parameter(get_parameter_reference(np, internal_np::weights_param), def_speeds);
    return SSEI::extrude_skeleton(pwh, weights, out, np);
  }
  else if(has_angles)
  {
    // not using "Weights" here is on purpose, I want to copy the range but the NP is ref only
    auto angles = choose_parameter(get_parameter_reference(np, internal_np::angles_param), def_speeds);
    SSEI::convert_angles<FT>(angles);
    return SSEI::extrude_skeleton(pwh, angles, out, np);
  }
  else // neither angles nor weights were provided
  {
    std::vector<std::vector<FT> > uniform_weights;
    uniform_weights.reserve(pwh.number_of_holes() + 1);

    uniform_weights.push_back(std::vector<FT>(pwh.outer_boundary().size(), FT(1)));
    for(const auto& hole : pwh.holes())
      uniform_weights.push_back(std::vector<FT>(hole.size(), FT(1)));

    return SSEI::extrude_skeleton(pwh, uniform_weights, out, np);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Polygon, // without holes
          typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
bool extrude_skeleton(const Polygon& p,
                      PolygonMesh& out,
                      const NamedParameters& np = parameters::default_values(),
                      std::enable_if_t<!CGAL_SS_i::has_Hole_const_iterator<Polygon>::value>* = nullptr)
{
  using Polygon_with_holes = CGAL::General_polygon_with_holes_2<Polygon>;

  return extrude_skeleton(Polygon_with_holes(p), out, np);
}

} // namespace CGAL

#endif // CGAL_SLS_EXTRUSION_EXTRUDE_SKELETON_H
