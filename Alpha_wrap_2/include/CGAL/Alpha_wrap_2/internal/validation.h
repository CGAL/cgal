// Copyright (c) 2025 GeometryFactory (France).
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
#ifndef CGAL_ALPHA_WRAP_2_TEST_ALPHA_WRAP_VALIDATION_H
#define CGAL_ALPHA_WRAP_2_TEST_ALPHA_WRAP_VALIDATION_H

#include <CGAL/license/Alpha_wrap_2.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Polygon_repair/repair.h>

#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

template <typename MultipolygonWithHoles>
bool is_empty(const MultipolygonWithHoles& wrap)
{
  using Polygon_with_holes_2 = typename MultipolygonWithHoles::Polygon_with_holes_2;

  if (wrap.polygons_with_holes().empty()) {
    return true;
  }

  // The real logic would be that it's empty if everything is empty,
  // but an empty pwh in the range is an error, so the logic here is
  // to mark it as empty as soon as one polygon with holes is empty.
  for (const Polygon_with_holes_2& pwh : wrap.polygons_with_holes()) {
    if(pwh.outer_boundary().size() == 0) {
      return true;
    }
  }

  return false;
}

template <typename MultipolygonWithHoles>
auto cdt_from_wrap(const MultipolygonWithHoles& wrap)
{
  using Polygon_with_holes_2 = typename MultipolygonWithHoles::Polygon_with_holes_2;
  using Polygon_2 = typename Polygon_with_holes_2::Polygon_2;
  using Point_2 = typename Polygon_2::Point_2;
  using K = typename CGAL::Kernel_traits<Point_2>::type;

  using Itag = CGAL::No_constraint_intersection_tag;
  using CDT = Constrained_triangulation_2<K, CGAL::Default, Itag>;
  using CDT_FH = typename CDT::Face_handle;

  CDT cdt;

  try
  {
    for (const Polygon_with_holes_2& pwh : wrap.polygons_with_holes()) {
      cdt.insert_constraint(pwh.outer_boundary().vertices_begin(),
                            pwh.outer_boundary().vertices_end(),
                            true /*close*/);

      for (const Polygon_2& hole : pwh.holes()) {
        cdt.insert_constraint(hole.vertices_begin(), hole.vertices_end(), true /*close*/);
      }
    }
  }
  catch(const typename CDT::Intersection_of_constraints_exception&)
  {
    CGAL_assertion_msg(false, "Intersections in CDT2 of wrap are not allowed");
    return CDT{};
  }

  return cdt;
}

// This is not a generic weakly simple multipolygon detection algorithm: we know that
// the multipolygon is a wrap and was constructed as the boundary of a set of inside cells.
// As such, it cannot twist, only touch, and it's enough to check that there are
// no intersections between edges which are not a subface of the edges.
template <typename MultipolygonWithHoles>
bool is_weakly_simple(const MultipolygonWithHoles& wrap)
{
  auto cdt = cdt_from_wrap(wrap);
  return cdt.number_of_vertices() != 0;
}

template <typename MultipolygonWithHoles>
bool is_simple(const MultipolygonWithHoles& wrap)
{
  auto cdt = cdt_from_wrap(wrap);

  using CDT = decltype(cdt);
  using Edge_circulator = typename CDT::Edge_circulator;

  for (auto vh : cdt.finite_vertex_handles()) {
    Edge_circulator ec = cdt.incident_edges(vh), done(ec);
    unsigned int cntr = 0;
    do {
      if(cdt.is_constrained(*ec)) {
        ++cntr;
      }
    } while (++ec != done);

    if (cntr != 2) {
#ifdef CGAL_AW2_DEBUG
      std::cerr << "Vertex " << cdt.point(vh) << " has " << cntr << " incident constraints" << std::endl;
#endif
      return false;
    }
  }

  return true;
}

// Edge length is bounded by twice the circumradius
template <typename MultipolygonWithHoles>
bool has_bounded_edge_length(const MultipolygonWithHoles& wrap,
                             const double alpha)
{
  using Polygon_with_holes_2 = typename MultipolygonWithHoles::Polygon_with_holes_2;
  using Polygon_2 = typename Polygon_with_holes_2::Polygon_2;
  using Edges = typename Polygon_2::Edges;
  using Segment_2 = typename Polygon_2::Segment_2;

  const auto sq_alpha_bound = 4 * square(alpha);

  for (const Polygon_with_holes_2& pwh : wrap.polygons_with_holes()) {
    for (const Segment_2& edge : pwh.outer_boundary().edges()) {
      if (CGAL::compare_squared_distance(edge.source(), edge.target(), sq_alpha_bound) == CGAL::LARGER) {
        return true;
      }
    }
    for (const Polygon_2& hole : pwh.holes()) {
      for (const Segment_2& edge : hole.edges()) {
        if (CGAL::compare_squared_distance(edge.source(), edge.target(), sq_alpha_bound) == CGAL::LARGER) {
          return true;
        }
      }
    }
  }

  return false;
}

// template <typename ConcurrencyTag = CGAL::Sequential_tag,
//           typename MultipolygonWithHoles, typename FT,
//           typename InputNamedParameters = parameters::Default_named_parameters,
//           typename OutputNamedParameters = parameters::Default_named_parameters>
// bool has_expected_Hausdorff_distance(const MultipolygonWithHoles& wrap,
//                                      const TriangleMesh& input,
//                                      const FT alpha, const FT offset,
//                                      const InputNamedParameters& in_np = parameters::default_values(),
//                                      const OutputNamedParameters& out_np = parameters::default_values())
// {
//   // @todo
//   return true;
// }

template <typename MultipolygonWithHoles,
          typename NamedParameters = parameters::Default_named_parameters>
bool is_valid_wrap(const MultipolygonWithHoles& wrap,
                   const bool check_manifoldness,
                   const NamedParameters& np = parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  if(Alpha_wraps_2::internal::is_empty(wrap))
  {
#ifdef CGAL_AW2_DEBUG
    std::cerr << "Error: empty wrap" << std::endl;
#endif
    return false;
  }

  if(!Polygon_repair::is_valid(wrap))
  {
#ifdef CGAL_AW2_DEBUG
    std::cerr << "Error: invalid wrap" << std::endl;
#endif
    return false;
  }

  if(check_manifoldness)
  {
    if(!Alpha_wraps_2::internal::is_simple(wrap))
    {
#ifdef CGAL_AW2_DEBUG
      std::cerr << "Error: Wrap is not simple" << std::endl;
#endif
    }
  }
  else
  {
    if(!Alpha_wraps_2::internal::is_weakly_simple(wrap))
    {
#ifdef CGAL_AW2_DEBUG
      std::cerr << "Error: Wrap is not weakly simple" << std::endl;
#endif
    }
  }

  return true;
}

template <typename Point, typename CDT, typename InDomain>
bool is_point_inside_marked_CDT(const Point& p,
                                const CDT& cdt,
                                const InDomain in_domain)
{
  using CDT_FH = typename CDT::Face_handle;
  using CDT_VH = typename CDT::Vertex_handle;
  using Locate_type = typename CDT::Locate_type;
  using Face_circulator = typename CDT::Face_circulator;

  CDT_FH fh;
  int li;
  Locate_type lt;
  fh = cdt.locate(p, lt, li);

  if(lt == CDT::VERTEX) {
    CDT_VH vh = fh->vertex(li);
    bool is_in = false;
    Face_circulator fc = cdt.incident_faces(vh), done(fc);
    do {
      if (get(in_domain, fc)) {
        is_in = true;
        break;
      }
    } while(++fc != done);
    if (!is_in) {
#ifdef CGAL_AW2_DEBUG
      std::cerr << "A polyline point [on vertex] is outside the wrap: " << p << std::endl;
#endif
      return false;
    }
  } else if (lt == CDT::EDGE) {
    if (!get(in_domain, fh) && !get(in_domain, fh->neighbor(li))) {
#ifdef CGAL_AW2_DEBUG
      std::cerr << "A polyline point [on edge] is outside the wrap: " << p << std::endl;
#endif
      return false;
    }
  } else if (lt == CDT::FACE) {
    if (!get(in_domain, fh)) {
#ifdef CGAL_AW2_DEBUG
      std::cerr << "A polyline point [on face] is outside the wrap: " << p << std::endl;
#endif
      return false;
    }
  } else {
#ifdef CGAL_AW2_DEBUG
    std::cerr << "A polyline point is outside of convex/affine hull?! " << p << std::endl;
#endif
    return false;
  }

  return true;
}

//////////////////////////////////////// POLYLINES /////////////////////////////////////////////////

template <typename MultipolygonWithHoles, typename PolylineRange,
          typename OutputNamedParameters = parameters::Default_named_parameters,
          typename InputNamedParameters = parameters::Default_named_parameters>
bool is_outer_wrap_of_polylines(const MultipolygonWithHoles& wrap,
                                const PolylineRange& polylines,
                                const OutputNamedParameters& out_np = parameters::default_values(),
                                const InputNamedParameters& in_np = parameters::default_values())
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using IPM = typename GetPointMap<PolylineRange, InputNamedParameters>::const_type;
  IPM in_pm = choose_parameter<IPM>(get_parameter(in_np, internal_np::point_map));

  using Polygon_with_holes_2 = typename MultipolygonWithHoles::Polygon_with_holes_2;
  using Polygon_2 = typename Polygon_with_holes_2::Polygon_2;
  using Segment_2 = typename Polygon_2::Segment_2;
  using Point_2 = typename Polygon_2::Point_2;
  using Box = CGAL::Box_intersection_d::Box_with_handle_d<double, 2, Segment_2*>;

  // Build boxed segments for wrap
  std::vector<Segment_2> wrap_segments;
  for(const Polygon_with_holes_2& pwh : wrap.polygons_with_holes()) {
    for(auto eit = pwh.outer_boundary().edges_begin(); eit != pwh.outer_boundary().edges_end(); ++eit)
      wrap_segments.push_back(*eit);
    for(const Polygon_2& hole : pwh.holes())
      for(auto eit = hole.edges_begin(); eit != hole.edges_end(); ++eit)
        wrap_segments.push_back(*eit);
  }
  std::vector<Box> wrap_boxes;
  wrap_boxes.reserve(wrap_segments.size());
  for(auto& seg : wrap_segments) {
    wrap_boxes.emplace_back(seg.bbox(), &seg);
  }

  // Build boxed segments for input polylines
  std::vector<Segment_2> input_segments;
  for(const auto& polyline : polylines) {
    if(polyline.size() < 2)
      continue;
    for(std::size_t i=1; i<polyline.size(); ++i)
      input_segments.emplace_back(polyline[i-1], polyline[i]);
  }
  std::vector<Box> input_boxes;
  input_boxes.reserve(input_segments.size());
  for(Segment_2& seg : input_segments) {
    input_boxes.emplace_back(seg.bbox(), &seg);
  }

  // Intersection callback
  struct FirstIntersection {};
  struct Overlap {
    void operator()(const Box& a, const Box& b) const {
      if(CGAL::do_intersect(*(a.handle()), *(b.handle())))
        throw FirstIntersection();
    }
  };

  try {
    CGAL::box_intersection_d(wrap_boxes.begin(), wrap_boxes.end(),
                             input_boxes.begin(), input_boxes.end(),
                             Overlap());
  } catch(const FirstIntersection&) {
#ifdef CGAL_AW2_DEBUG
    std::cerr << "Polyline segment intersects wrap boundary" << std::endl;
#endif
    return false;
  }

  // Mark domain
  auto cdt = cdt_from_wrap(wrap);

  using CDT = decltype(cdt);
  using CDT_FH = typename CDT::Face_handle;

  std::unordered_map<CDT_FH, bool> in_domain_map;
  boost::associative_property_map<std::unordered_map<CDT_FH, bool> > in_domain(in_domain_map);
  CGAL::mark_domain_in_triangulation(cdt, in_domain);

  for(const auto& polyline : polylines) {
    for(const auto& p : polyline) {
      if(!is_point_inside_marked_CDT(p, cdt, in_domain))
        return false;
    }
  }

  return true;
}

template <typename MultipolygonWithHoles, typename PolylineRange,
          typename OutputNamedParameters = parameters::Default_named_parameters,
          typename InputNamedParameters = parameters::Default_named_parameters>
bool is_valid_wrap_of_polylines(const MultipolygonWithHoles& wrap,
                                const PolylineRange& polylines,
                                const OutputNamedParameters& out_np = parameters::default_values(),
                                const InputNamedParameters& in_np = parameters::default_values())
{
  if(!is_valid_wrap(wrap, out_np))
    return false;

  if(!is_outer_wrap_of_polylines(wrap, polylines, out_np, in_np))
    return false;

  return true;
}

//////////////////////////////////////// POINT SET /////////////////////////////////////////////////

template <typename MultipolygonWithHoles, typename PointRange,
          typename OutputNamedParameters = parameters::Default_named_parameters,
          typename InputNamedParameters = parameters::Default_named_parameters>
bool is_outer_wrap_of_point_set(const MultipolygonWithHoles& wrap,
                                const PointRange& points,
                                const OutputNamedParameters& out_np = parameters::default_values(),
                                const InputNamedParameters& in_np = parameters::default_values())
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using IPM = typename GetPointMap<PointRange, InputNamedParameters>::const_type;
  using K = typename Kernel_traits<typename boost::property_traits<IPM>::value_type>::Kernel;

  IPM in_pm = choose_parameter<IPM>(get_parameter(in_np, internal_np::point_map));

  auto cdt = cdt_from_wrap(wrap);

  using CDT = decltype(cdt);
  using CDT_VH = typename CDT::Vertex_handle;
  using CDT_FH = typename CDT::Face_handle;
  using Locate_type = typename CDT::Locate_type;
  using Face_circulator = typename CDT::Face_circulator;

  std::unordered_map<CDT_FH, bool> in_domain_map;
  boost::associative_property_map<std::unordered_map<CDT_FH, bool> > in_domain(in_domain_map);
  CGAL::mark_domain_in_triangulation(cdt, in_domain);

  for(const auto& p : points)
    if(!is_point_inside_marked_CDT(get(in_pm, p), cdt, in_domain))
      return false;

  return true;
}

template <typename MultipolygonWithHoles, typename PointRange,
          typename OutputNamedParameters = parameters::Default_named_parameters,
          typename InputNamedParameters = parameters::Default_named_parameters>
bool is_valid_wrap_of_point_set(const MultipolygonWithHoles& wrap,
                                const PointRange& points,
                                const OutputNamedParameters& out_np = parameters::default_values(),
                                const InputNamedParameters& in_np = parameters::default_values())
{
  if(!is_valid_wrap(wrap, out_np))
    return false;

  if(!is_outer_wrap_of_point_set(wrap, points, out_np, in_np))
    return false;

  return true;
}

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_TEST_ALPHA_WRAP_VALIDATION_H
