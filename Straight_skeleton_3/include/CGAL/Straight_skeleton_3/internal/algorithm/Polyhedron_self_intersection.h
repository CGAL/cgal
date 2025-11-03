// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   algo/3d/SelfIntersection.h
 * author Gernot Walzl
 * date   2012-07-18
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_SELF_INTERSECTION_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_SELF_INTERSECTION_H

#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>

#include <CGAL/enum.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Random.h>

#include <limits>
#include <list>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename GeomTraits>
class Self_intersection
{
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;
  using Segment_3 = typename GeomTraits::Segment_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  using Ray_3 = typename GeomTraits::Ray_3;
  using Line_3 = typename GeomTraits::Line_3;
  using Plane_3 = typename GeomTraits::Plane_3;

private:
  using Polyhedron = HDS::Polyhedron<GeomTraits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

public:
  // check if edges share at least one vertex
  static bool do_edges_share_a_vertex(const EdgeSPtr& edge1,
                                      const EdgeSPtr& edge2,
                                      bool handle_degree_1_as_ray)
  {
    CGAL_SS3_DEBUG_SPTR(edge1);
    CGAL_SS3_DEBUG_SPTR(edge2);

    // if handle_degree_1_as_ray is true, then that extremitiy is not considered
    VertexSPtr v1_src = edge1->source();
    if (handle_degree_1_as_ray && v1_src->degree() == 1) {
      v1_src = nullptr;
    }
    VertexSPtr v1_tgt = edge1->target();
    if (handle_degree_1_as_ray && v1_tgt->degree() == 1) {
      v1_tgt = nullptr;
    }
    VertexSPtr v2_src = edge2->source();
    if (handle_degree_1_as_ray && v2_src->degree() == 1) {
      v2_src = nullptr;
    }
    VertexSPtr v2_tgt = edge2->target();
    if (handle_degree_1_as_ray && v2_tgt->degree() == 1) {
      v2_tgt = nullptr;
    }

    return (v1_src == v2_src || v1_src == v2_tgt ||
            v1_tgt == v2_src || v1_tgt == v2_tgt);
  }

  // whether the edges intersect in their interior, with edges possibly being rays
  template <typename ProjectionTraits>
  static bool do_facet_edges_intersect(const FacetSPtr&,
                                       const EdgeSPtr& edge1,
                                       const EdgeSPtr& edge2,
                                       bool handle_degree_1_as_ray, // @todo always true?...
                                       const ProjectionTraits& traits)
  {
    CGAL_SS3_DEBUG_SPTR(edge1);
    CGAL_SS3_DEBUG_SPTR(edge2);
    CGAL_precondition(edge1 != edge2);
    CGAL_precondition(edge1->source() != edge1->target() &&
                      edge2->source() != edge2->target());

    // edges should never be a full line
    CGAL_precondition(edge1->source()->degree() != 1 || edge1->target()->degree() != 1);
    CGAL_precondition(edge2->source()->degree() != 1 || edge2->target()->degree() != 1);

    // reject any degeneracy as a self-intersection
    if (edge1->source()->point() == edge1->target()->point() ||
        edge2->source()->point() == edge2->target()->point()) {
      return false;
    }

    // we could have something like an overlapping segment and ray, but in that case,
    // we will see the intersection from another pair
    if (do_edges_share_a_vertex(edge1, edge2, handle_degree_1_as_ray)) {
      return false;
    }

    auto is_ray = [](const EdgeSPtr& edge) -> bool
    {
      return (edge->source()->degree() == 1 || edge->target()->degree() == 1);
    };

    auto treat_o_o = [&](const EdgeSPtr&, const auto& oa,
                         const EdgeSPtr&, const auto& ob) -> bool
    {
      return traits.do_intersect_2_object()(oa, ob);
    };

    auto treat_seg_seg = [&](const EdgeSPtr& a, const EdgeSPtr& b) -> bool
    {
      CGAL_precondition(!is_ray(a) && !is_ray(b));
      const Segment_3 sa = { a->source()->point(),
                             a->target()->point() };
      const Segment_3 sb = { b->source()->point(),
                             b->target()->point() };
      return treat_o_o(a, sa, b, sb);
    };

    auto treat_seg_ray = [&](const EdgeSPtr& a, const EdgeSPtr& b) -> bool
    {
      CGAL_precondition(!is_ray(a) && is_ray(b));
      const Segment_3 s = { a->source()->point(),
                            a->target()->point() };
      const Ray_3 r = (b->source()->degree() == 1) ? Ray_3 { b->target()->point(),
                                                             b->source()->point() }
                                                   : Ray_3 { b->source()->point(),
                                                             b->target()->point() };
      return treat_o_o(a, s, b, r);
    };

    auto treat_ray_ray = [&](const EdgeSPtr& a, const EdgeSPtr& b) -> bool
    {
      CGAL_precondition(is_ray(a) && is_ray(b));
      const Ray_3 ra = (a->source()->degree() == 1) ? Ray_3 { a->target()->point(),
                                                              a->source()->point() }
                                                    : Ray_3 { a->source()->point(),
                                                              a->target()->point() };
      const Ray_3 rb = (b->source()->degree() == 1) ? Ray_3 { b->target()->point(),
                                                              b->source()->point() }
                                                    : Ray_3 { b->source()->point(),
                                                              b->target()->point() };
      return treat_o_o(a, ra, b, rb);
    };

    if (handle_degree_1_as_ray) {
      if (is_ray(edge1)) {
        if (is_ray(edge2)) {
          return treat_ray_ray(edge1, edge2);
        } else {
          return treat_seg_ray(edge2, edge1);
        }
      } else if (is_ray(edge2)) {
        return treat_seg_ray(edge1, edge2);
      }
    }

    return treat_seg_seg(edge1, edge2);
  }

  template <typename ProjectionTraits>
  static bool is_self_intersecting_facet(const FacetSPtr& facet,
                                         const ProjectionTraits& traits)
  {
    CGAL_SS3_DEBUG_SPTR(facet);

    bool result = false;

    // We can't use is_simple_polygon_2 + projection traits because some edges are rays
    typename std::list<EdgeSPtr>::iterator it_e1 = facet->edges().begin();
    while (it_e1 != facet->edges().end()) {
      EdgeSPtr edge1 = *it_e1++;
      typename std::list<EdgeSPtr>::iterator it_e2 = it_e1;
      while (it_e2 != facet->edges().end()) {
        EdgeSPtr edge2 = *it_e2++;
        if (do_facet_edges_intersect(facet, edge1, edge2, true, traits)) {
          CGAL_SS3_ALGO_TRACE_V(32, "edges intersect within the facet:");
          CGAL_SS3_ALGO_TRACE_V(32, facet->to_string());
          CGAL_SS3_ALGO_TRACE_V(32, edge1->to_string());
          CGAL_SS3_ALGO_TRACE_V(32, edge2->to_string());
#ifdef CGAL_SS3_EXIT_ASAP
          return true;
#else
          result = true;
#endif
          break;
        }
      }
      if (result) {
          break;
      }
    }
    return result;
  }

  static bool is_self_intersecting_facet(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    const Plane_3& pl = facet->get_plane();
    const Vector_3 normal = pl.orthogonal_vector();
    if (is_zero(normal.z())) {
      using Traits_2 = CGAL::Projection_traits_xz_3<GeomTraits>;
      Traits_2 traits;
      return is_self_intersecting_facet(facet, traits);
    } else {
      using Traits_2 = CGAL::Projection_traits_xy_3<GeomTraits>;
      Traits_2 traits;
      return is_self_intersecting_facet(facet, traits);
    }
  }

  static unsigned int has_self_intersecting_facets(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    unsigned int result = 0;
    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (is_self_intersecting_facet(facet)) {
#ifdef CGAL_SS3_EXIT_ASAP
        return 1;
#else
        ++result;
#endif
      }
    }
    return result;
  }

  // @todo construction galore...
  static bool is_inside_using_ray_shooting(const Point_3& point,
                                           const FacetSPtr& facet,
                                           const bool handle_degree_1_as_ray)
  {
    CGAL_SS3_DEBUG_SPTR(facet);

    CGAL_SS3_ALGO_TRACE_V(64, "is_inside() using ray shooting (" << point << ", F" << facet->id() << ")");

    const Plane_3& pl = facet->get_plane();
    const Vector_3 normal = pl.orthogonal_vector();

    // shoot random rays till something is hit
    // essential to this: we know the facet does not self-intersect
    CGAL::Random rng(0);
    CGAL::Random_points_on_sphere_3<Point_3> random_point_on_sphere(1, rng);

    // some initial rays that we know are in the plane and will hit something
    std::vector<Point_3> candidate_ray_targets;

    for (const EdgeSPtr& edge : facet->edges()) {
      const Point_3& p_src = edge->source()->point();
      const Point_3& p_tgt = edge->target()->point();

      CGAL_SS3_ALGO_TRACE_V(64, " p_src = " << p_src);
      CGAL_SS3_ALGO_TRACE_V(64, " p_tgt = " << p_tgt);
      CGAL_assertion(p_src != p_tgt);

      // this only stands for EPECK and flat faces
      CGAL_assertion(pl.has_on(point));
      CGAL_assertion(pl.has_on(p_src));
      CGAL_assertion(pl.has_on(p_tgt));
      CGAL_assertion(CGAL::scalar_product(Vector_3(p_src, p_tgt), normal) == 0);

      if (handle_degree_1_as_ray) {
        CGAL_assertion(edge->source() != edge->target());
        CGAL_precondition(edge->source()->degree() != 1 || edge->target()->degree() != 1);

        VertexSPtr r_src = nullptr;
        VertexSPtr r_tgt = nullptr;
        if (edge->source()->degree() == 1) {
          r_src = edge->target();
          r_tgt = edge->source();
        } else if (edge->target()->degree() == 1) {
          r_src = edge->source();
          r_tgt = edge->target();
        }

        if (r_src && r_tgt) {
          Ray_3 r { r_src->point(), r_tgt->point() };
          if (r.has_on(point)) {
            // intersection if it's on the ray except if its the source
            return (point != r_src->point());
          }
        } else {
          Segment_3 s { p_src, p_tgt };
          if (s.has_on(point)) {
            // intersection if it's on the ray except if its the source or target
            return (point != p_src && point != p_tgt);
          }
        }
      } else {
        Segment_3 s { p_src, p_tgt };
        if (s.has_on(point)) {
          // intersection if it's on the ray except if its the source or target
          return (point != p_src && point != p_tgt);
        }
      }

      candidate_ray_targets.push_back(CGAL::midpoint(p_src, p_tgt));

      CGAL_SS3_ALGO_TRACE_V(64, "new potential ray target: " << candidate_ray_targets.back() << " for " << edge->to_string());
    }

    for (;;) {
      Ray_3 shooting_ray;
      if (!candidate_ray_targets.empty()) {
        shooting_ray = Ray_3(point, candidate_ray_targets.back());
        candidate_ray_targets.pop_back();
      } else {
        Point_3 rnd_p = *random_point_on_sphere++;
        Point_3 target_p = point + Vector_3(CGAL::ORIGIN, rnd_p);
        Point_3 proj_p = pl.projection(target_p);
        if (proj_p == point) {
          continue;
        }
        shooting_ray = Ray_3(point, proj_p);
      }

      CGAL_SS3_ALGO_TRACE_V(64, "shooting_ray = " << shooting_ray.point(0) << " " << shooting_ray.point(1));

      CGAL_assertion(shooting_ray.point(0) == point);
      CGAL_assertion(pl.has_on(shooting_ray.point(0)));
      CGAL_assertion(pl.has_on(shooting_ray.point(1)));

      // normally we shouldn't need any of the random targets
      if (shooting_ray.is_degenerate()) {
        continue;
      }

      FT sq_dist_to_closest = std::numeric_limits<double>::max();
      EdgeSPtr closest_edge;

      auto treat_edge = [&](const EdgeSPtr& target_edge, const auto& edge_geometry) -> void
      {
        CGAL_assertion(!edge_geometry.is_degenerate());

        CGAL_SS3_ALGO_TRACE_V(64, "Treat " << target_edge->to_string());

        CGAL::Object obj = CGAL::intersection(shooting_ray, edge_geometry);
        CGAL_assertion_code(using EG = CGAL::cpp20::remove_cvref_t<decltype(edge_geometry)>);
        CGAL_assertion_code(const EG* eg = CGAL::object_cast<EG>(&obj);)
        CGAL_assertion(!bool(eg));

        if (const Point_3 *ipoint = CGAL::object_cast<Point_3>(&obj)) {
          FT sqd = CGAL::squared_distance(point, *ipoint);
          CGAL_SS3_ALGO_TRACE_V(64, "  intersects @ SQ_tgt: " << sqd);
          CGAL::Comparison_result res = CGAL::compare(sqd, sq_dist_to_closest);
          if (res == CGAL::SMALLER) {
            sq_dist_to_closest = sqd;
            closest_edge = target_edge;
            CGAL_SS3_ALGO_TRACE_V(64, "  new closest");
          }
        }
      };

      // @speed this is a brute force approach, but apart from debugging, we only use
      // self-intersection detection while splitting or handling edge events,
      // where facets have few edges
      for (const EdgeSPtr& edge : facet->edges()) {
        CGAL_SS3_ALGO_TRACE_V(64, "consider edge: " << edge->to_string());
        VertexSPtr v_src = edge->src(facet);
        VertexSPtr v_tgt = edge->tgt(facet);
        const Point_3& p_src = v_src->point();
        const Point_3& p_tgt = v_tgt->point();

        if (CGAL::collinear(p_src, p_tgt, point)) {
          // we have already checked that the point is not on an edge
          // while collecting ray targets
          continue;
        }

        if (handle_degree_1_as_ray) {
          if (v_src->degree() == 1) {
            if (v_tgt->degree() == 1) {
              treat_edge(edge, Line_3(p_src, p_tgt));
            } else {
              treat_edge(edge, Ray_3(p_tgt, p_src));
            }
          } else if (v_tgt->degree() == 1) {
            treat_edge(edge, Ray_3(p_src, p_tgt));
          } else {
            treat_edge(edge, Segment_3(p_src, p_tgt));
          }
        } else {
          treat_edge(edge, Segment_3(p_src, p_tgt));
        }
      }

      CGAL_SS3_ALGO_TRACE_CODE(if (closest_edge))
      CGAL_SS3_ALGO_TRACE_V(64, "closest_edge = " << closest_edge->to_string());

      // Being in the cone where the distance is the same to multiple edges
      // makes the orientation test not usable
      if (closest_edge == EdgeSPtr()) {
        continue; // try another ray, we will hit something eventually!
      }

      const Point_3& p_src = closest_edge->src(facet)->point();
      const Point_3& p_tgt = closest_edge->tgt(facet)->point();

      CGAL_SS3_ALGO_TRACE_V(64, "p_src = " << p_src);
      CGAL_SS3_ALGO_TRACE_V(64, "p_src + normal = " << p_src + normal);
      CGAL_SS3_ALGO_TRACE_V(64, "p_tgt = " << p_tgt);
      CGAL_SS3_ALGO_TRACE_V(64, "point = " << point);

      CGAL_assertion(!CGAL::collinear(p_src, p_src + normal, p_tgt));
      CGAL_assertion(CGAL::scalar_product(Vector_3(p_src, p_src + normal), Vector_3(p_src, p_tgt)) == 0);

      CGAL::Orientation o = CGAL::orientation(p_src, p_src + normal, p_tgt, point);
      CGAL_SS3_ALGO_TRACE_V(64, "Orientation = " << o);

      return (o != CGAL::NEGATIVE);
    }

    CGAL_unreachable();
    return false;
  }

private:
  // returns -1 if point is left of segment <low, high>, 0 if its on the segment
  // and 1 if it is to the right
  // precondition: low.y < point.y < high.y
  template <class Point, class Orientation_2, class CompareX_2>
  static int which_side_in_slab(const Point& point,
                                const Point& low,
                                const Point& high,
                                Orientation_2& orientation_2,
                                CompareX_2& compare_x_2)
  {
    // first we try to decide on x coordinate values alone
    // This is an optimization (whether this is really faster for
    // a homogeneous kernel is not clear, as comparisons can be expensive.
    CGAL::Comparison_result low_x_comp_res = compare_x_2(point, low);
    CGAL::Comparison_result high_x_comp_res = compare_x_2(point, high);
    if (low_x_comp_res == CGAL::SMALLER) {
      if (high_x_comp_res == CGAL::SMALLER)
        return -1;
    } else {
      switch (high_x_comp_res) {
        case CGAL::LARGER: return 1;
        case CGAL::SMALLER: break;
        case  CGAL::EQUAL: return (low_x_comp_res ==  CGAL::EQUAL) ? 0 : 1;
      }
    }
    switch (orientation_2(low, point, high)) {
      case CGAL::LEFT_TURN: return 1;
      case CGAL::RIGHT_TURN: return -1;
      default: return 0;
    }
  }

  template <typename ProjectionTraits>
  static CGAL::Bounded_side bounded_side(const Point_3& point,
                                         const FacetSPtr& facet,
                                         const ProjectionTraits& traits)
  {
    bool is_inside = false;

    // Iterate over all edges, treating each as a segment in the projected plane
    for (const EdgeSPtr& edge : facet->edges()) {
      const Point_3& p_src = edge->src(facet)->point();
      const Point_3& p_tgt = edge->tgt(facet)->point();

      // Ray-shooting logic: check if the edge crosses the horizontal ray from point
      typename ProjectionTraits::Compare_y_2 compare_y_2 = traits.compare_y_2_object();
      typename ProjectionTraits::Compare_x_2 compare_x_2 = traits.compare_x_2_object();
      typename ProjectionTraits::Orientation_2 orientation_2 = traits.orientation_2_object();

      CGAL::Comparison_result src_y = compare_y_2(p_src, point);
      CGAL::Comparison_result tgt_y = compare_y_2(p_tgt, point);

      switch (src_y) {
        case CGAL::SMALLER:
          switch (tgt_y) {
            case CGAL::SMALLER:
              break;
            case  CGAL::EQUAL:
              switch (compare_x_2(point, p_tgt)) {
                case CGAL::SMALLER: is_inside = !is_inside; break;
                case  CGAL::EQUAL:   return CGAL::ON_BOUNDARY;
                case CGAL::LARGER:  break;
              }
              break;
            case CGAL::LARGER:
              switch (which_side_in_slab(point, p_src, p_tgt, orientation_2, compare_x_2)) {
                case -1: is_inside = !is_inside; break;
                case  0: return CGAL::ON_BOUNDARY;
              }
              break;
          }
          break;
        case  CGAL::EQUAL:
          switch (tgt_y) {
            case CGAL::SMALLER:
              switch (compare_x_2(point, p_src)) {
                case CGAL::SMALLER: is_inside = !is_inside; break;
                case  CGAL::EQUAL:   return CGAL::ON_BOUNDARY;
                case CGAL::LARGER:  break;
              }
              break;
            case  CGAL::EQUAL:
              switch (compare_x_2(point, p_src)) {
                case CGAL::SMALLER:
                  if (compare_x_2(point, p_tgt) != CGAL::SMALLER)
                      return CGAL::ON_BOUNDARY;
                  break;
                case  CGAL::EQUAL: return CGAL::ON_BOUNDARY;
                case CGAL::LARGER:
                  if (compare_x_2(point, p_tgt) != CGAL::LARGER)
                    return CGAL::ON_BOUNDARY;
                  break;
              }
              break;
            case CGAL::LARGER:
              if (compare_x_2(point, p_src) ==  CGAL::EQUAL) {
                return CGAL::ON_BOUNDARY;
              }
              break;
          }
          break;
        case CGAL::LARGER:
          switch (tgt_y) {
            case CGAL::SMALLER:
              switch (which_side_in_slab(point, p_tgt, p_src, orientation_2, compare_x_2)) {
                case -1: is_inside = !is_inside; break;
                case  0: return CGAL::ON_BOUNDARY;
              }
              break;
            case  CGAL::EQUAL:
              if (compare_x_2(point, p_tgt) ==  CGAL::EQUAL) {
                return CGAL::ON_BOUNDARY;
              }
              break;
            case CGAL::LARGER:
              break;
          }
          break;
      }
    }

    return is_inside ? CGAL::ON_BOUNDED_SIDE : CGAL::ON_UNBOUNDED_SIDE;
  }

public:
  static bool is_inside_using_ray_shooting_V2(const Point_3& point,
                                              const FacetSPtr& facet)
  {
    CGAL_SS3_ALGO_TRACE_V(32, "is_inside() using ray shooting V2 (" << point << ", F" << facet->id() << ")");

    const Plane_3& pl = facet->get_plane();
    const Vector_3 normal = pl.orthogonal_vector();
    if (is_zero(normal.z())) {
      using Traits_2 = CGAL::Projection_traits_xz_3<GeomTraits>;
      Traits_2 traits;
      return (bounded_side(point, facet, traits) != CGAL::ON_UNBOUNDED_SIDE);
    } else {
      using Traits_2 = CGAL::Projection_traits_xy_3<GeomTraits>;
      Traits_2 traits;
      return (bounded_side(point, facet, traits) != CGAL::ON_UNBOUNDED_SIDE);
    }
  }

  static bool is_edge_inside_facet(const FacetSPtr& facet,
                                const EdgeSPtr& edge,
                                bool handle_degree_1_as_ray)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_SS3_DEBUG_SPTR(edge);

    CGAL_SS3_ALGO_TRACE_V(32, "\n> is edge inside facet");
    CGAL_SS3_ALGO_TRACE_V(32, "  " << facet->to_string());
    CGAL_SS3_ALGO_TRACE_V(32, "  " << edge->to_string());

    VertexSPtr e_src = edge->source();
    VertexSPtr e_tgt = edge->target();

    CGAL_precondition(e_src != e_tgt);

    // edges should never be a full line
    CGAL_precondition(e_src->degree() != 1 || e_tgt->degree() != 1);

    if (edge->get_facet_L() == facet || edge->get_facet_R() == facet) {
      return false;
    }

    const Plane_3& facet_pl = facet->get_plane();

    // Start with the case of the edge living in the same plane as the facet (unlikely to happen)
    // @speed use the projection traits here too
    bool coplanarity = (facet->contains_vertex(e_src) || facet_pl.has_on(e_src->point())) &&
                       (facet->contains_vertex(e_tgt) || facet_pl.has_on(e_tgt->point()));
    if (coplanarity) {
      auto test_fo_o_coplanarity = [&](const EdgeSPtr& /*fe*/, const auto& fo,
                                       const EdgeSPtr& /*e*/, const auto& o) -> bool
      {
        // - fe and e do not share a vertex
        // - if the intersection is 1-dimensional, it's a real intersection because 'edge'
        // is not incident to the facet
        return CGAL::do_intersect(fo, o);
      };

      auto test_fe_o_coplanarity = [&](const EdgeSPtr& fe, const EdgeSPtr& e, const auto& o) -> bool
      {
        // we could have something like an overlapping segment and ray, but in that case,
        // we will see the intersection from another pair
        if (do_edges_share_a_vertex(fe, e, handle_degree_1_as_ray)) {
          return false;
        }

        // distinguish if the _facet edge_ is a ray or a segment
        if (handle_degree_1_as_ray) {
          VertexSPtr r_src = nullptr;
          VertexSPtr r_tgt = nullptr;
          if (fe->source()->degree() == 1) {
            r_src = fe->target();
            r_tgt = fe->source();
          } else if (fe->target()->degree() == 1) {
            r_src = fe->source();
            r_tgt = fe->target();
          }

          if (r_src && r_tgt) {
            Ray_3 r { r_src->point(), r_tgt->point() };
            return test_fo_o_coplanarity(fe, r, e, o);
          }
        }

        Segment_3 s { fe->source()->point(),
                      fe->target()->point() };
        return test_fo_o_coplanarity(fe, s, e, o);
      };

      auto test_edges = [&](const EdgeSPtr& e, const auto& o) -> bool
      {
        for (const EdgeSPtr& facet_edge : facet->edges()) {
          if (test_fe_o_coplanarity(facet_edge, e, o)) {
            return true;
          }
        }

        return false;
      };

      // distinguish if the _querying edge_ is a ray or a segment
      if (handle_degree_1_as_ray) {
        VertexSPtr r_src = nullptr;
        VertexSPtr r_tgt = nullptr;
        if (e_src->degree() == 1) {
          r_src = e_tgt;
          r_tgt = e_src;
        } else if (e_tgt->degree() == 1) {
          r_src = e_src;
          r_tgt = e_tgt;
        }

        if (r_src && r_tgt) {
          Ray_3 r { r_src->point(), r_tgt->point() };
          if (test_edges(edge, r)) {
            return true;
          }
        } else {
          Segment_3 s { e_src->point(), e_tgt->point() };
          if (test_edges(edge, s)) {
            return true;
          }
        }
      } else {
        Segment_3 s { e_src->point(), e_tgt->point() };
        if (test_edges(edge, s)) {
          return true;
        }
      }

      // if there is no intersection between the edge and any facet edge, then an extremity
      // of the edge is sufficient to determine where we are
      return is_inside_using_ray_shooting(e_src->point(), facet, handle_degree_1_as_ray);
    }

    // Here we know that the edge does not live in the plane of the facet

    auto test_o_facet = [&](const EdgeSPtr& e, const auto& o) -> bool
    {
      using O = CGAL::cpp20::remove_cvref_t<decltype(o)>;

      CGAL_SS3_ALGO_TRACE_V(64, "test_o_facet(" << e->id() << " " << typeid(O).name() << ")");

      CGAL::Object obj = CGAL::intersection(facet_pl, o);
      if (!obj) {
        return false; // no intersection
      } else if (const Point_3 *ipoint = CGAL::object_cast<Point_3>(&obj)) {
        CGAL_SS3_ALGO_TRACE_V(64, "intersection point at " << *ipoint);

        // intersection is a point, so there is an intersection if the point is not an extremity
        // of the edge, and if it is inside the facet
        if constexpr (std::is_same_v<O, Segment_3>) {
          if (facet->contains_vertex(e->source()) ||
              facet->contains_vertex(e->target())) {
            return false;
          } else {
              // @speed adapt v2 to work with rays?
            return is_inside_using_ray_shooting(*ipoint, facet, handle_degree_1_as_ray);
          }
        } else if constexpr (std::is_same_v<O, Ray_3>) {
          VertexSPtr r_src = (e->source()->degree() == 1) ? e->target()
                                                          : e->source();
          if (facet->contains_vertex(r_src)) {
            return false;
          } else {
            return is_inside_using_ray_shooting(*ipoint, facet, handle_degree_1_as_ray);
          }
        }
      } else {
        // Here, the intersection is 1-dimensional. This shouldn't be possible
        // because the edge is not coplanar with the facet.
        CGAL_assertion(false);
        return true;
      }
    };

    if (handle_degree_1_as_ray) {
      VertexSPtr r_src = nullptr;
      VertexSPtr r_tgt = nullptr;
      if (e_src->degree() == 1) {
        r_src = e_tgt;
        r_tgt = e_src;
      } else if (e_tgt->degree() == 1) {
        r_src = e_src;
        r_tgt = e_tgt;
      }

      if (r_src && r_tgt) {
        Ray_3 r { r_src->point(), r_tgt->point() };
        return test_o_facet(edge, r);
      }
    }

    Segment_3 s { e_src->point(), e_tgt->point() };
    return test_o_facet(edge, s);
  }

  static bool has_self_intersecting_surface(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    CGAL_SS3_ALGO_TRACE_V(16, "has_self_intersecting_surface()");

    // @speed O(nf*nfe*nfe) algorithm
    if (Self_intersection::has_self_intersecting_facets(polyhedron)) {
      return true;
    }

    // @speed O(nf*nfe*ne) algorithm
    for (const FacetSPtr& facet : polyhedron->facets()) {
      for (const EdgeSPtr& edge : polyhedron->edges()) {
        if (is_edge_inside_facet(facet, edge, true /*handle_degree_1_as_ray*/)) {
          CGAL_SS3_ALGO_TRACE_V(16, "\nPolyhedron has no self-intersecting facets, but the surface is self-intersecting!");
          CGAL_SS3_ALGO_TRACE_V(16, facet->to_string());
          CGAL_SS3_ALGO_TRACE_V(16, edge->to_string());
          return true;
        }
      }
    }

    return false;
  }
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_SELF_INTERSECTION_H */
