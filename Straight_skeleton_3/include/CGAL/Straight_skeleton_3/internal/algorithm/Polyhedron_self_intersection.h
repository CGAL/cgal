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
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>

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

template <typename Traits>
class SelfIntersection
{
  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Segment_3 = typename Traits::Segment_3;
  using Vector_3 = typename Traits::Vector_3;
  using Ray_3 = typename Traits::Ray_3;
  using Line_3 = typename Traits::Line_3;
  using Plane_3 = typename Traits::Plane_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

private:
  using KernelFactory = kernel::KernelFactory<Traits>;

public:
  // check if edges share at least one vertex
  static bool doEdgesShareAVertex(const EdgeSPtr& edge1,
                                  const EdgeSPtr& edge2,
                                  bool handle_degree_1_as_ray)
  {
    CGAL_SS3_DEBUG_SPTR(edge1);
    CGAL_SS3_DEBUG_SPTR(edge2);

    // if handle_degree_1_as_ray is true, then that extremitiy is not considered
    VertexSPtr v1_src = edge1->getVertexSrc();
    if (handle_degree_1_as_ray && v1_src->degree() == 1) {
      v1_src = nullptr;
    }
    VertexSPtr v1_dst = edge1->getVertexDst();
    if (handle_degree_1_as_ray && v1_dst->degree() == 1) {
      v1_dst = nullptr;
    }
    VertexSPtr v2_src = edge2->getVertexSrc();
    if (handle_degree_1_as_ray && v2_src->degree() == 1) {
      v2_src = nullptr;
    }
    VertexSPtr v2_dst = edge2->getVertexDst();
    if (handle_degree_1_as_ray && v2_dst->degree() == 1) {
      v2_dst = nullptr;
    }

    return (v1_src == v2_src || v1_src == v2_dst ||
            v1_dst == v2_src || v1_dst == v2_dst);
  }

  // whether the edges intersect in their interior, with edges possibly being rays
  template <typename ProjectionTraits>
  static bool doEdgesIntersect(const FacetSPtr&,
                               const EdgeSPtr& edge1,
                               const EdgeSPtr& edge2,
                               bool handle_degree_1_as_ray, // @todo always true?...
                               const ProjectionTraits& traits)
  {
    CGAL_SS3_DEBUG_SPTR(edge1);
    CGAL_SS3_DEBUG_SPTR(edge2);
    CGAL_precondition(edge1 != edge2);
    CGAL_precondition(edge1->getVertexSrc() != edge1->getVertexDst() &&
                      edge2->getVertexSrc() != edge2->getVertexDst());

    // edges should never be a full line
    CGAL_precondition(edge1->getVertexSrc()->degree() != 1 || edge1->getVertexDst()->degree() != 1);
    CGAL_precondition(edge2->getVertexSrc()->degree() != 1 || edge2->getVertexDst()->degree() != 1);

    // reject any degeneracy as a self-intersection
    if (*(edge1->getVertexSrc()->getPoint()) == *(edge1->getVertexDst()->getPoint()) ||
        *(edge2->getVertexSrc()->getPoint()) == *(edge2->getVertexDst()->getPoint())) {
      return false;
    }

    // we could have something like an overlapping segment and ray, but in that case,
    // we will see the intersection from another pair
    if (doEdgesShareAVertex(edge1, edge2, handle_degree_1_as_ray)) {
      return false;
    }

    auto isRay = [](const EdgeSPtr& edge) -> bool
    {
      return (edge->getVertexSrc()->degree() == 1 || edge->getVertexDst()->degree() == 1);
    };

    auto treat_o_o = [&](const EdgeSPtr&, const auto& oa,
                         const EdgeSPtr&, const auto& ob) -> bool
    {
      return traits.do_intersect_2_object()(oa, ob);
    };

    auto treat_seg_seg = [&](const EdgeSPtr& a, const EdgeSPtr& b) -> bool
    {
      CGAL_precondition(!isRay(a) && !isRay(b));
      const Segment_3 sa = { *(a->getVertexSrc()->getPoint()),
                             *(a->getVertexDst()->getPoint()) };
      const Segment_3 sb = { *(b->getVertexSrc()->getPoint()),
                             *(b->getVertexDst()->getPoint()) };
      return treat_o_o(a, sa, b, sb);
    };

    auto treat_seg_ray = [&](const EdgeSPtr& a, const EdgeSPtr& b) -> bool
    {
      CGAL_precondition(!isRay(a) && isRay(b));
      const Segment_3 s = { *(a->getVertexSrc()->getPoint()),
                            *(a->getVertexDst()->getPoint()) };
      const Ray_3 r = (b->getVertexSrc()->degree() == 1) ? Ray_3 { *b->getVertexDst()->getPoint(),
                                                                   *b->getVertexSrc()->getPoint() }
                                                         : Ray_3 { *b->getVertexSrc()->getPoint(),
                                                                   *b->getVertexDst()->getPoint() };
      return treat_o_o(a, s, b, r);
    };

    auto treat_ray_ray = [&](const EdgeSPtr& a, const EdgeSPtr& b) -> bool
    {
      CGAL_precondition(isRay(a) && isRay(b));
      const Ray_3 ra = (a->getVertexSrc()->degree() == 1) ? Ray_3 { *a->getVertexDst()->getPoint(),
                                                                    *a->getVertexSrc()->getPoint() }
                                                          : Ray_3 { *a->getVertexSrc()->getPoint(),
                                                                    *a->getVertexDst()->getPoint() };
      const Ray_3 rb = (b->getVertexSrc()->degree() == 1) ? Ray_3 { *b->getVertexDst()->getPoint(),
                                                                    *b->getVertexSrc()->getPoint() }
                                                          : Ray_3 { *b->getVertexSrc()->getPoint(),
                                                                    *b->getVertexDst()->getPoint() };
      return treat_o_o(a, ra, b, rb);
    };

    if (handle_degree_1_as_ray) {
      if (isRay(edge1)) {
        if (isRay(edge2)) {
          return treat_ray_ray(edge1, edge2);
        } else {
          return treat_seg_ray(edge2, edge1);
        }
      } else if (isRay(edge2)) {
        return treat_seg_ray(edge1, edge2);
      }
    }

    return treat_seg_seg(edge1, edge2);
  }

  template <typename ProjectionTraits>
  static bool isSelfIntersectingFacet(const FacetSPtr& facet,
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
        if (doEdgesIntersect(facet, edge1, edge2, true, traits)) {
          CGAL_SS3_ALGO_TRACE_V(32, "edges intersect within the facet:");
          CGAL_SS3_ALGO_TRACE_V(32, facet->toString());
          CGAL_SS3_ALGO_TRACE_V(32, edge1->toString());
          CGAL_SS3_ALGO_TRACE_V(32, edge2->toString());
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

  static bool isSelfIntersectingFacet(const FacetSPtr& facet)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    Plane3SPtr pl = facet->plane();
    Vector3SPtr normal = KernelFactory::createVector3(pl);
    if (is_zero(normal->z())) {
      typedef CGAL::Projection_traits_xz_3<Traits> Traits_2;
      Traits_2 traits;
      return isSelfIntersectingFacet(facet, traits);
    } else {
      typedef CGAL::Projection_traits_xy_3<Traits> Traits_2;
      Traits_2 traits;
      return isSelfIntersectingFacet(facet, traits);
    }
  }

  static unsigned int hasSelfIntersectingFacets(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    unsigned int result = 0;
    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (isSelfIntersectingFacet(facet)) {
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
  static bool isInsideWithRayShooting(const Point_3& point,
                                      const FacetSPtr& facet,
                                      const bool handle_degree_1_as_ray)
  {
    CGAL_SS3_DEBUG_SPTR(facet);

    CGAL_SS3_ALGO_TRACE_V(64, "isInsideWithRayShooting(" << point << ", F" << facet->getID() << ")");

    Plane3SPtr pl = facet->plane();
    Vector3SPtr normal = KernelFactory::createVector3(pl);

    // shoot random rays till something is hit
    // essential to this: we know the facet does not self-intersect
    CGAL::Random rng(0);
    CGAL::Random_points_on_sphere_3<Point_3> random_point_on_sphere(1, rng);

    // some initial rays that we know are in the plane and will hit something
    std::vector<Point_3> candidate_ray_targets;

    for (const EdgeSPtr& edge : facet->edges()) {
      Point3SPtr p_src = edge->getVertexSrc()->getPoint();
      Point3SPtr p_dst = edge->getVertexDst()->getPoint();

      CGAL_SS3_ALGO_TRACE_V(64, " p_src = " << *p_src);
      CGAL_SS3_ALGO_TRACE_V(64, " p_dst = " << *p_dst);
      CGAL_assertion(*p_src != *p_dst);

      // this only stands for EPECK and flat faces
      CGAL_assertion(pl->has_on(point));
      CGAL_assertion(pl->has_on(*p_src));
      CGAL_assertion(pl->has_on(*p_dst));
      CGAL_assertion(CGAL::scalar_product(Vector_3(*p_src, *p_dst), *normal) == 0);

      if (handle_degree_1_as_ray) {
        CGAL_assertion(edge->getVertexSrc() != edge->getVertexDst());
        CGAL_precondition(edge->getVertexSrc()->degree() != 1 || edge->getVertexDst()->degree() != 1);

        VertexSPtr r_src = nullptr;
        VertexSPtr r_dst = nullptr;
        if (edge->getVertexSrc()->degree() == 1) {
          r_src = edge->getVertexDst();
          r_dst = edge->getVertexSrc();
        } else if (edge->getVertexDst()->degree() == 1) {
          r_src = edge->getVertexSrc();
          r_dst = edge->getVertexDst();
        }

        if (r_src && r_dst) {
          Ray_3 r { *r_src->getPoint(), *r_dst->getPoint() };
          if (r.has_on(point)) {
            // intersection if it's on the ray except if its the source
            return (point != *(r_src->getPoint()));
          }
        } else {
          Segment_3 s { *p_src, *p_dst };
          if (s.has_on(point)) {
            // intersection if it's on the ray except if its the source or target
            return (point != *p_src && point != *p_dst);
          }
        }
      } else {
        Segment_3 s { *p_src, *p_dst };
        if (s.has_on(point)) {
          // intersection if it's on the ray except if its the source or target
          return (point != *p_src && point != *p_dst);
        }
      }

      candidate_ray_targets.push_back(CGAL::midpoint(*p_src, *p_dst));

      CGAL_SS3_ALGO_TRACE_V(64, "new potential ray target: " << candidate_ray_targets.back() << " for " << edge->toString());
    }

    for (;;) {
      Ray_3 shooting_ray;
      if (!candidate_ray_targets.empty()) {
        shooting_ray = Ray_3(point, candidate_ray_targets.back());
        candidate_ray_targets.pop_back();
      } else {
        Point_3 rnd_p = *random_point_on_sphere++;
        Point_3 target_p = point + Vector_3(CGAL::ORIGIN, rnd_p);
        Point_3 proj_p = pl->projection(target_p);
        if (proj_p == point) {
          continue;
        }
        shooting_ray = Ray_3(point, proj_p);
      }

      CGAL_SS3_ALGO_TRACE_V(64, "shooting_ray = " << shooting_ray.point(0) << " " << shooting_ray.point(1));

      CGAL_assertion(shooting_ray.point(0) == point);
      CGAL_assertion(pl->has_on(shooting_ray.point(0)));
      CGAL_assertion(pl->has_on(shooting_ray.point(1)));

      // normally we shouldn't need any of the random targets
      if (shooting_ray.is_degenerate()) {
        continue;
      }

      FT sq_dist_to_closest = std::numeric_limits<double>::max();
      EdgeSPtr closest_edge;

      auto treat_edge = [&](const EdgeSPtr& target_edge, const auto& edge_geometry) -> void
      {
        CGAL_assertion(!edge_geometry.is_degenerate());

        CGAL_SS3_ALGO_TRACE_V(64, "Treat " << target_edge->toString());

        CGAL::Object obj = CGAL::intersection(shooting_ray, edge_geometry);
        CGAL_assertion_code(using EG = CGAL::cpp20::remove_cvref_t<decltype(edge_geometry)>);
        CGAL_assertion_code(const EG* eg = CGAL::object_cast<EG>(&obj);)
        CGAL_assertion(!bool(eg));

        if (const Point_3 *ipoint = CGAL::object_cast<Point_3>(&obj)) {
          FT sqd = CGAL::squared_distance(point, *ipoint);
          CGAL_SS3_ALGO_TRACE_V(64, "  intersects @ SQ_dst: " << sqd);
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
        CGAL_SS3_ALGO_TRACE_V(64, "consider edge: " << edge->toString());
        VertexSPtr v_src = edge->src(facet);
        VertexSPtr v_dst = edge->dst(facet);
        Point3SPtr p_src = v_src->getPoint();
        Point3SPtr p_dst = v_dst->getPoint();

        if (CGAL::collinear(*p_src, *p_dst, point)) {
          // we have already checked that the point is not on an edge
          // while collecting ray targets
          continue;
        }

        if (handle_degree_1_as_ray) {
          if (v_src->degree() == 1) {
            if (v_dst->degree() == 1) {
              treat_edge(edge, Line_3(*p_src, *p_dst));
            } else {
              treat_edge(edge, Ray_3(*p_dst, *p_src));
            }
          } else if (v_dst->degree() == 1) {
            treat_edge(edge, Ray_3(*p_src, *p_dst));
          } else {
            treat_edge(edge, Segment_3(*p_src, *p_dst));
          }
        } else {
          treat_edge(edge, Segment_3(*p_src, *p_dst));
        }
      }

      CGAL_SS3_ALGO_TRACE_CODE(if (closest_edge))
      CGAL_SS3_ALGO_TRACE_V(64, "closest_edge = " << closest_edge->toString());

      // Being in the cone where the distance is the same to multiple edges
      // makes the orientation test not usable
      if (closest_edge == EdgeSPtr()) {
        continue; // try another ray, we will hit something eventually!
      }

      Point3SPtr p_src = closest_edge->src(facet)->getPoint();
      Point3SPtr p_dst = closest_edge->dst(facet)->getPoint();

      CGAL_SS3_ALGO_TRACE_V(64, "p_src = " << *p_src);
      CGAL_SS3_ALGO_TRACE_V(64, "p_src + normal = " << *p_src + *normal);
      CGAL_SS3_ALGO_TRACE_V(64, "p_dst = " << *p_dst);
      CGAL_SS3_ALGO_TRACE_V(64, "point = " << point);

      CGAL_assertion(!CGAL::collinear(*p_src, *p_src + *normal, *p_dst));
      CGAL_assertion(CGAL::scalar_product(Vector_3(*p_src, *p_src + *normal), Vector_3(*p_src, *p_dst)) == 0);

      CGAL::Orientation o = CGAL::orientation(*p_src, *p_src + *normal, *p_dst, point);
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
  static CGAL::Bounded_side boundedSide(const Point3SPtr& point,
                                        const FacetSPtr& facet,
                                        const ProjectionTraits& traits)
  {
    bool is_inside = false;

    // Iterate over all edges, treating each as a segment in the projected plane
    for (const EdgeSPtr& edge : facet->edges()) {
      Point3SPtr p_src = edge->src(facet)->getPoint();
      Point3SPtr p_dst = edge->dst(facet)->getPoint();

      // Ray-shooting logic: check if the edge crosses the horizontal ray from point
      typename ProjectionTraits::Compare_y_2 compare_y_2 = traits.compare_y_2_object();
      typename ProjectionTraits::Compare_x_2 compare_x_2 = traits.compare_x_2_object();
      typename ProjectionTraits::Orientation_2 orientation_2 = traits.orientation_2_object();

      CGAL::Comparison_result src_y = compare_y_2(*p_src, *point);
      CGAL::Comparison_result dst_y = compare_y_2(*p_dst, *point);

      switch (src_y) {
        case CGAL::SMALLER:
          switch (dst_y) {
            case CGAL::SMALLER:
              break;
            case  CGAL::EQUAL:
              switch (compare_x_2(*point, *p_dst)) {
                case CGAL::SMALLER: is_inside = !is_inside; break;
                case  CGAL::EQUAL:   return CGAL::ON_BOUNDARY;
                case CGAL::LARGER:  break;
              }
              break;
            case CGAL::LARGER:
              switch (which_side_in_slab(*point, *p_src, *p_dst, orientation_2, compare_x_2)) {
                case -1: is_inside = !is_inside; break;
                case  0: return CGAL::ON_BOUNDARY;
              }
              break;
          }
          break;
        case  CGAL::EQUAL:
          switch (dst_y) {
            case CGAL::SMALLER:
              switch (compare_x_2(*point, *p_src)) {
                case CGAL::SMALLER: is_inside = !is_inside; break;
                case  CGAL::EQUAL:   return CGAL::ON_BOUNDARY;
                case CGAL::LARGER:  break;
              }
              break;
            case  CGAL::EQUAL:
              switch (compare_x_2(*point, *p_src)) {
                case CGAL::SMALLER:
                  if (compare_x_2(*point, *p_dst) != CGAL::SMALLER)
                      return CGAL::ON_BOUNDARY;
                  break;
                case  CGAL::EQUAL: return CGAL::ON_BOUNDARY;
                case CGAL::LARGER:
                  if (compare_x_2(*point, *p_dst) != CGAL::LARGER)
                    return CGAL::ON_BOUNDARY;
                  break;
              }
              break;
            case CGAL::LARGER:
              if (compare_x_2(*point, *p_src) ==  CGAL::EQUAL) {
                return CGAL::ON_BOUNDARY;
              }
              break;
          }
          break;
        case CGAL::LARGER:
          switch (dst_y) {
            case CGAL::SMALLER:
              switch (which_side_in_slab(*point, *p_dst, *p_src, orientation_2, compare_x_2)) {
                case -1: is_inside = !is_inside; break;
                case  0: return CGAL::ON_BOUNDARY;
              }
              break;
            case  CGAL::EQUAL:
              if (compare_x_2(*point, *p_dst) ==  CGAL::EQUAL) {
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
  static bool isInsideWithRayShootingV2(const Point3SPtr& point,
                                        const FacetSPtr& facet)
  {
    CGAL_SS3_ALGO_TRACE_V(32, "isInsideWithRayShootingV2(" << point << ", F" << facet->getID() << ")");

    Plane3SPtr pl = facet->plane();
    Vector3SPtr normal = KernelFactory::createVector3(pl);
    if (is_zero(normal->z())) {
      typedef CGAL::Projection_traits_xz_3<Traits> Traits_2;
      Traits_2 traits;
      return (boundedSide(point, facet, traits) != CGAL::ON_UNBOUNDED_SIDE);
    } else {
      typedef CGAL::Projection_traits_xy_3<Traits> Traits_2;
      Traits_2 traits;
      return (boundedSide(point, facet, traits) != CGAL::ON_UNBOUNDED_SIDE);
    }
  }

  static bool isEdgeInsideFacet(const FacetSPtr& facet,
                                const EdgeSPtr& edge,
                                bool handle_degree_1_as_ray)
  {
    CGAL_SS3_DEBUG_SPTR(facet);
    CGAL_SS3_DEBUG_SPTR(edge);

    CGAL_SS3_ALGO_TRACE_V(32, "\n> isEdgeInsideFacet()");
    CGAL_SS3_ALGO_TRACE_V(32, "  " << facet->toString());
    CGAL_SS3_ALGO_TRACE_V(32, "  " << edge->toString());

    VertexSPtr e_src = edge->getVertexSrc();
    VertexSPtr e_dst = edge->getVertexDst();

    CGAL_precondition(e_src != e_dst);

    // edges should never be a full line
    CGAL_precondition(e_src->degree() != 1 || e_dst->degree() != 1);

    if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
      return false;
    }

    Plane3SPtr facet_pl = facet->plane();

    // Start with the case of the edge living in the same plane as the facet (unlikely to happen)
    // @speed use the projection traits here too
    bool coplanarity = (facet->containsVertex(e_src) || facet_pl->has_on(*e_src->getPoint())) &&
                       (facet->containsVertex(e_dst) || facet_pl->has_on(*e_dst->getPoint()));
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
        if (doEdgesShareAVertex(fe, e, handle_degree_1_as_ray)) {
          return false;
        }

        // distinguish if the _facet edge_ is a ray or a segment
        if (handle_degree_1_as_ray) {
          VertexSPtr r_src = nullptr;
          VertexSPtr r_dst = nullptr;
          if (fe->getVertexSrc()->degree() == 1) {
            r_src = fe->getVertexDst();
            r_dst = fe->getVertexSrc();
          } else if (fe->getVertexDst()->degree() == 1) {
            r_src = fe->getVertexSrc();
            r_dst = fe->getVertexDst();
          }

          if (r_src && r_dst) {
            Ray_3 r { *r_src->getPoint(), *r_dst->getPoint() };
            return test_fo_o_coplanarity(fe, r, e, o);
          }
        }

        Segment_3 s { *fe->getVertexSrc()->getPoint(),
                      *fe->getVertexDst()->getPoint() };
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
        VertexSPtr r_dst = nullptr;
        if (e_src->degree() == 1) {
          r_src = e_dst;
          r_dst = e_src;
        } else if (e_dst->degree() == 1) {
          r_src = e_src;
          r_dst = e_dst;
        }

        if (r_src && r_dst) {
          Ray_3 r { *r_src->getPoint(), *r_dst->getPoint() };
          if (test_edges(edge, r)) {
            return true;
          }
        } else {
          Segment_3 s { *e_src->getPoint(), *e_dst->getPoint() };
          if (test_edges(edge, s)) {
            return true;
          }
        }
      } else {
        Segment_3 s { *e_src->getPoint(), *e_dst->getPoint() };
        if (test_edges(edge, s)) {
          return true;
        }
      }

      // if there is no intersection between the edge and any facet edge, then an extremity
      // of the edge is sufficient to determine where we are
      return isInsideWithRayShooting(*(e_src->getPoint()), facet, handle_degree_1_as_ray);
    }

    // Here we know that the edge does not live in the plane of the facet

    auto test_o_facet = [&](const EdgeSPtr& e, const auto& o) -> bool
    {
      using O = CGAL::cpp20::remove_cvref_t<decltype(o)>;

      CGAL_SS3_ALGO_TRACE_V(64, "test_o_facet(" << e->getID() << " " << typeid(O).name() << ")");

      CGAL::Object obj = CGAL::intersection(*facet_pl, o);
      if (!obj) {
        return false; // no intersection
      } else if (const Point_3 *ipoint = CGAL::object_cast<Point_3>(&obj)) {
        CGAL_SS3_ALGO_TRACE_V(64, "intersection point at " << *ipoint);

        // intersection is a point, so there is an intersection if the point is not an extremity
        // of the edge, and if it is inside the facet
        if constexpr (std::is_same_v<O, Segment_3>) {
          if (facet->containsVertex(e->getVertexSrc()) ||
              facet->containsVertex(e->getVertexDst())) {
            return false;
          } else {
              // @speed adapt v2 to work with rays?
            return isInsideWithRayShooting(*ipoint, facet, handle_degree_1_as_ray);
          }
        } else if constexpr (std::is_same_v<O, Ray_3>) {
          VertexSPtr r_src = (e->getVertexSrc()->degree() == 1) ? e->getVertexDst()
                                                                : e->getVertexSrc();
          if (facet->containsVertex(r_src)) {
            return false;
          } else {
            return isInsideWithRayShooting(*ipoint, facet, handle_degree_1_as_ray);
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
      VertexSPtr r_dst = nullptr;
      if (e_src->degree() == 1) {
        r_src = e_dst;
        r_dst = e_src;
      } else if (e_dst->degree() == 1) {
        r_src = e_src;
        r_dst = e_dst;
      }

      if (r_src && r_dst) {
        Ray_3 r { *r_src->getPoint(), *r_dst->getPoint() };
        return test_o_facet(edge, r);
      }
    }

    Segment_3 s { *e_src->getPoint(), *e_dst->getPoint() };
    return test_o_facet(edge, s);
  }

  // @fixme self-intersection on boundaries does not seem well defined
  static bool hasSelfIntersectingSurface(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    CGAL_SS3_ALGO_TRACE_V(16, "hasSelfIntersectingSurface()");

    // @speed O(nf*nfe*nfe) algorithm
    if (SelfIntersection::hasSelfIntersectingFacets(polyhedron)) {
      return true;
    }

    // @speed O(nf*nfe*ne) algorithm
    for (const FacetSPtr& facet : polyhedron->facets()) {
      for (const EdgeSPtr& edge : polyhedron->edges()) {
        if (isEdgeInsideFacet(facet, edge, true /*handle_degree_1_as_ray*/)) {
          CGAL_SS3_ALGO_TRACE_V(16, "\nPolyhedron has no self-intersecting facets, but the surface is self-intersecting!");
          CGAL_SS3_ALGO_TRACE_V(16, facet->toString());
          CGAL_SS3_ALGO_TRACE_V(16, edge->toString());
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
