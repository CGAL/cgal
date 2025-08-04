// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   algo/3d/SelfIntersection.cpp
 * @author Gernot Walzl
 * @date   2012-07-18
 */

#include "algo/3d/SelfIntersection.h"

#include "data/3d/KernelFactory.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "algo/3d/KernelWrapper.h"

#include <CGAL/point_generators_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_xz_3.h>

#include <limits>
#include <list>

#define CGAL_SS3_EXIT_ASAP

namespace algo { namespace _3d {

// check if edges share a vertex (could be more, but we only use this functor
// when there is a single point of intersection, so it cannot be more)
bool SelfIntersection::doEdgesShareAVertex(EdgeSPtr edge1,
                                           EdgeSPtr edge2,
                                           bool handle_degree_1_as_ray)
{
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
};

// whether the edges intersect in their interior, with edges possibly being rays
bool SelfIntersection::doEdgesIntersect(FacetSPtr facet,
                                        EdgeSPtr edge1,
                                        EdgeSPtr edge2,
                                        bool handle_degree_1_as_ray)
{
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

    auto isRay = [](EdgeSPtr edge) -> bool
    {
        return (edge->getVertexSrc()->degree() == 1 || edge->getVertexDst()->degree() == 1);
    };

    auto treat_o_o = [&](EdgeSPtr a, const auto& oa,
                         EdgeSPtr b, const auto& ob) -> bool
    {
        CGAL::Object obj = CGAL::intersection(oa, ob);
        if (!obj) {
            return false;
        } else if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
            // intersection is a point; there is an intersection unless that point is a common extremity
            // (note that it needs to be a shared vertex, and not just the same point
            return !doEdgesShareAVertex(a, b, handle_degree_1_as_ray);
        } else {
            return true; // intersection is a segment
        }
    };

    auto treat_seg_seg = [&](EdgeSPtr a, EdgeSPtr b) -> bool
    {
        CGAL_precondition(!isRay(a) && !isRay(b));
        const Segment3 sa = { *(a->getVertexSrc()->getPoint()),
                              *(a->getVertexDst()->getPoint()) };
        const Segment3 sb = { *(b->getVertexSrc()->getPoint()),
                              *(b->getVertexDst()->getPoint()) };
        return treat_o_o(a, sa, b, sb);
    };

    auto treat_seg_ray = [&](EdgeSPtr a, EdgeSPtr b) -> bool
    {
        CGAL_precondition(!isRay(a) && isRay(b));
        const Segment3 s = { *(a->getVertexSrc()->getPoint()),
                             *(a->getVertexDst()->getPoint()) };
        const Ray3 r = (b->getVertexSrc()->degree() == 1) ? Ray3 { *b->getVertexDst()->getPoint(),
                                                                   *b->getVertexSrc()->getPoint() }
                                                          : Ray3 { *b->getVertexSrc()->getPoint(),
                                                                   *b->getVertexDst()->getPoint() };
        return treat_o_o(a, s, b, r);
    };

    auto treat_ray_ray = [&](EdgeSPtr a, EdgeSPtr b) -> bool
    {
        CGAL_precondition(isRay(a) && isRay(b));
        const Ray3 ra = (a->getVertexSrc()->degree() == 1) ? Ray3 { *a->getVertexDst()->getPoint(),
                                                                    *a->getVertexSrc()->getPoint() }
                                                           : Ray3 { *a->getVertexSrc()->getPoint(),
                                                                    *a->getVertexDst()->getPoint() };
        const Ray3 rb = (b->getVertexSrc()->degree() == 1) ? Ray3 { *b->getVertexDst()->getPoint(),
                                                                    *b->getVertexSrc()->getPoint() }
                                                           : Ray3 { *b->getVertexSrc()->getPoint(),
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

// @speed currently has square complexity, it could be a sweep like in "is_simple_polygon_2"
bool SelfIntersection::isSelfIntersectingFacet(FacetSPtr facet) {
    bool result = false;

    // We can't use is_simple_polygon_2 + projection traits because:
    // - some edges are not segments but ray
    // - some edges can be degenerate
    std::list<EdgeSPtr>::iterator it_e1 = facet->edges().begin();
    while (it_e1 != facet->edges().end()) {
        EdgeSPtr edge1 = *it_e1++;
        std::list<EdgeSPtr>::iterator it_e2 = it_e1;
        while (it_e2 != facet->edges().end()) {
            EdgeSPtr edge2 = *it_e2++;
            if (doEdgesIntersect(facet, edge1, edge2, true)) {
                CGAL_SS3_ALGO_TRACE("edges intersect within the face:");
                CGAL_SS3_ALGO_TRACE(facet->toString());
                CGAL_SS3_ALGO_TRACE(edge1->toString());
                CGAL_SS3_ALGO_TRACE(edge2->toString());
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

unsigned int SelfIntersection::hasSelfIntersectingFacets(PolyhedronSPtr polyhedron) {
    unsigned int result = 0;
    for (FacetSPtr facet : polyhedron->facets()) {
        if (isSelfIntersectingFacet(facet)) {
#ifdef CGAL_SS3_EXIT_ASAP
            return 1;
#else
            result++;
#endif
        }
    }
    return result;
}

// @todo construction galore...
bool SelfIntersection::isInsideWithRayShooting(const Point3& point,
                                               FacetSPtr facet,
                                               const bool handle_degree_1_as_ray)
{
    CGAL_SS3_ALGO_TRACE("isInsideWithRayShooting(" << point << ", F" << facet->getID() << ")");

    Plane3SPtr pl = facet->plane();
    Vector3SPtr normal = KernelFactory::createVector3(pl);

    // shoot random rays till something is hit
    // essential to this: we know the facet does not self-intersect
    CGAL::Random rng(0);
    CGAL::Random_points_on_sphere_3<Point3> random_point_on_sphere(1, rng);

    // some initial rays that we know are in the plane and will hit something
    std::vector<Point3> candidate_ray_targets;

    for (EdgeSPtr edge : facet->edges()) {
        Point3SPtr p_src = edge->getVertexSrc()->getPoint();
        Point3SPtr p_dst = edge->getVertexDst()->getPoint();

        CGAL_SS3_ALGO_TRACE(" p_src = " << *p_src);
        CGAL_SS3_ALGO_TRACE(" p_dst = " << *p_dst);
        CGAL_assertion(*p_src != *p_dst);

        // this only stands for EPECK and flat faces
        CGAL_assertion(pl->has_on(point));
        CGAL_assertion(pl->has_on(*p_src));
        CGAL_assertion(pl->has_on(*p_dst));
        CGAL_assertion(CGAL::scalar_product(Vector3(*p_src, *p_dst), *normal) == 0);

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
                Ray3 r { *r_src->getPoint(), *r_dst->getPoint() };
                if (r.has_on(point)) {
                    // intersection if it's on the ray except if its the source
                    return (point != *(r_src->getPoint()));
                }
            } else {
                Segment3 s { *p_src, *p_dst };
                if (s.has_on(point)) {
                    // intersection if it's on the ray except if its the source or target
                    return (point != *p_src && point != *p_dst);
                }
            }
        } else {
            Segment3 s { *p_src, *p_dst };
            if (s.has_on(point)) {
                // intersection if it's on the ray except if its the source or target
                return (point != *p_src && point != *p_dst);
            }
        }

        candidate_ray_targets.push_back(CGAL::midpoint(*p_src, *p_dst));

        CGAL_SS3_ALGO_TRACE("new potential ray target: " << candidate_ray_targets.back() << " for " << edge->toString());
    }

    for(;;) {
        Ray3 shooting_ray;
        if (!candidate_ray_targets.empty()) {
            shooting_ray = Ray3(point, candidate_ray_targets.back());
            candidate_ray_targets.pop_back();
        } else {
            Point3 rnd_p = *random_point_on_sphere++;
            Point3 target_p = point + Vector3(CGAL::ORIGIN, rnd_p);
            Point3 proj_p = pl->projection(target_p);
            if(proj_p == point) {
                continue;
            }
            shooting_ray = Ray3(point, proj_p);
        }

        CGAL_SS3_ALGO_TRACE("shooting_ray = " << shooting_ray.point(0) << " " << shooting_ray.point(1));

        CGAL_assertion(shooting_ray.point(0) == point);
        CGAL_assertion(pl->has_on(shooting_ray.point(0)));
        CGAL_assertion(pl->has_on(shooting_ray.point(1)));

        // normally we shouldn't need any of the random targets
        if (shooting_ray.is_degenerate()) {
            continue;
        }

        CGAL::FT sq_dist_to_closest = std::numeric_limits<double>::max();
        EdgeSPtr closest_edge;

        auto treat_edge = [&](EdgeSPtr target_edge, const auto& edge_geometry) -> void
        {
            CGAL_assertion(!edge_geometry.is_degenerate());

            CGAL_SS3_ALGO_TRACE("Treat " << target_edge->toString());

            CGAL::Object obj = CGAL::intersection(shooting_ray, edge_geometry);
            CGAL_assertion_code(using EG = CGAL::cpp20::remove_cvref_t<decltype(edge_geometry)>);
            CGAL_assertion_code(const EG* eg = CGAL::object_cast<EG>(&obj);)
            CGAL_assertion(!bool(eg));

            if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
                CGAL::FT sqd = CGAL::squared_distance(point, *ipoint);
                CGAL_SS3_ALGO_TRACE("  intersects @ SQ_dst: " << sqd);
                CGAL::Comparison_result res = CGAL::compare(sqd, sq_dist_to_closest);
                if (res == CGAL::SMALLER) {
                    sq_dist_to_closest = sqd;
                    closest_edge = target_edge;
                    CGAL_SS3_ALGO_TRACE("  new closest");
                }
            }
        };

        // @speed this is a brute force approach, but apart from debugging, we only use
        // self-intersection detection while splitting or handling edge events,
        // where facets have few edges
        for (EdgeSPtr edge : facet->edges()) {
            CGAL_SS3_ALGO_TRACE("consider edge: " << edge->toString());
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
                        treat_edge(edge, Line3(*p_src, *p_dst));
                    } else {
                        treat_edge(edge, Ray3(*p_dst, *p_src));
                    }
                } else if (v_dst->degree() == 1) {
                    treat_edge(edge, Ray3(*p_src, *p_dst));
                } else {
                    treat_edge(edge, Segment3(*p_src, *p_dst));
                }
            } else {
                treat_edge(edge, Segment3(*p_src, *p_dst));
            }
        }

        CGAL_SS3_ALGO_TRACE_CODE(if (closest_edge))
        CGAL_SS3_ALGO_TRACE("closest_edge = " << closest_edge->toString());

        // Being in the cone where the distance is the same to multiple edges
        // makes the orientation test not usable
        if (closest_edge == EdgeSPtr()) {
            continue; // try another ray, we will hit something eventually!
        }

        Point3SPtr p_src = closest_edge->src(facet)->getPoint();
        Point3SPtr p_dst = closest_edge->dst(facet)->getPoint();

        CGAL_SS3_ALGO_TRACE("p_src = " << *p_src);
        CGAL_SS3_ALGO_TRACE("p_src + normal = " << *p_src + *normal);
        CGAL_SS3_ALGO_TRACE("p_dst = " << *p_dst);
        CGAL_SS3_ALGO_TRACE("point = " << point);

        CGAL_assertion(!CGAL::collinear(*p_src, *p_src + *normal, *p_dst));
        CGAL_assertion(CGAL::scalar_product(Vector3(*p_src, *p_src + *normal), Vector3(*p_src, *p_dst)) == 0);

        CGAL::Orientation o = CGAL::orientation(*p_src, *p_src + *normal, *p_dst, point);
        CGAL_SS3_ALGO_TRACE("Orientation = " << o);

        return (o != CGAL::NEGATIVE);
    }

    CGAL_unreachable();
    return false;
}

// returns -1 if point is left of segment <low, high>, 0 if its on the segment
// and 1 if it is to the right
// precondition: low.y < point.y < high.y
template <class Point, class Orientation_2, class CompareX_2>
int which_side_in_slab(const Point& point,
                       const Point& low, const Point& high,
                       Orientation_2& orientation_2, CompareX_2& compare_x_2)
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
CGAL::Bounded_side boundedSide(Point3SPtr point, FacetSPtr facet, const ProjectionTraits& traits)
{
    bool is_inside = false;

    // Iterate over all edges, treating each as a segment in the projected plane
    for (EdgeSPtr edge : facet->edges()) {
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

bool SelfIntersection::isInsideWithRayShootingV2(Point3SPtr point,
                                                 FacetSPtr facet)
{
    CGAL_SS3_ALGO_TRACE("isInsideWithRayShootingV2(" << point << ", F" << facet->getID() << ")");
    Plane3SPtr pl = facet->plane();
    Vector3SPtr normal = KernelFactory::createVector3(pl);
    if (is_zero(normal->z())) {
        typedef CGAL::Projection_traits_xz_3<CGAL::K> Traits_2;
        Traits_2 traits;
        return (boundedSide(point, facet, traits) != CGAL::ON_UNBOUNDED_SIDE);
    } else {
        typedef CGAL::Projection_traits_xy_3<CGAL::K> Traits_2;
        Traits_2 traits;
        return (boundedSide(point, facet, traits) != CGAL::ON_UNBOUNDED_SIDE);
    }
}

bool SelfIntersection::isEdgeInsideFacet(FacetSPtr facet,
                                         EdgeSPtr edge,
                                         bool handle_degree_1_as_ray)
{
    CGAL_SS3_ALGO_TRACE("\n> isEdgeInsideFacet()");
    CGAL_SS3_ALGO_TRACE("  " << facet->toString());
    CGAL_SS3_ALGO_TRACE("  " << edge->toString());

    VertexSPtr e_src = edge->getVertexSrc();
    VertexSPtr e_dst = edge->getVertexDst();

    CGAL_precondition(e_src != e_dst);

    // edges should never be a full line
    CGAL_precondition(e_src->degree() != 1 || e_dst->degree() != 1);

    if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
        return false;
    }

    Plane3SPtr facet_pl = facet->plane();

    // Start with the case of the edge living in the same plane as the facet
    bool coplanarity = (facet->containsVertex(e_src) || facet_pl->has_on(*e_src->getPoint())) &&
                       (facet->containsVertex(e_dst) || facet_pl->has_on(*e_dst->getPoint()));

    if (coplanarity) {
        auto test_fo_o_coplanarity = [&](EdgeSPtr fe, const auto& fo,
                                         EdgeSPtr e, const auto& o) -> bool
        {
            CGAL::Object obj = CGAL::intersection(fo, o);
            if (!obj) {
                return false;
            } else if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
                // intersection is a point; there is an intersection unless that point is a common extremity
                // (note that it needs to be a shared vertex, and not just the same point
                return !doEdgesShareAVertex(fe, e, handle_degree_1_as_ray);
            } else {
                // intersection is 1-dimensional, so it's a real intersection because 'edge'
                // is not incident to the facet
                return true;
            }
        };

        auto test_fe_o_coplanarity = [&](EdgeSPtr fe, EdgeSPtr e, const auto& o) -> bool
        {
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
                    Ray3 r { *r_src->getPoint(), *r_dst->getPoint() };
                    return test_fo_o_coplanarity(fe, r, e, o);
                }
            }

            Segment3 s { *fe->getVertexSrc()->getPoint(),
                         *fe->getVertexDst()->getPoint() };
            return test_fo_o_coplanarity(fe, s, e, o);
        };

        auto test_edges = [&](EdgeSPtr e, const auto& o) -> bool
        {
            for (EdgeSPtr facet_edge : facet->edges()) {
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
                Ray3 r { *r_src->getPoint(), *r_dst->getPoint() };
                if (test_edges(edge, r)) {
                    return true;
                }
            } else {
                Segment3 s { *e_src->getPoint(), *e_dst->getPoint() };
                if (test_edges(edge, s)) {
                    return true;
                }
            }
        } else {
            Segment3 s { *e_src->getPoint(), *e_dst->getPoint() };
            if (test_edges(edge, s)) {
                return true;
            }
        }

        // if there is no intersection between the edge and any facet edge, then an extremity
        // of the edge is sufficient to determine where we are
        return isInsideWithRayShooting(*(e_src->getPoint()), facet, handle_degree_1_as_ray);
    }

    // Now we know that the edge does not live in the plane of the facet

    auto test_o_facet = [&](EdgeSPtr e, const auto& o) -> bool
    {
        using O = CGAL::cpp20::remove_cvref_t<decltype(o)>;

        CGAL_SS3_ALGO_TRACE("test_o_facet(" << e->getID() << " " << typeid(O).name() << ")");

        CGAL::Object obj = CGAL::intersection(*facet_pl, o);
        if (!obj) {
            return false; // no intersection
        } else if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
            CGAL_SS3_ALGO_TRACE("intersection point at " << *ipoint);

            // intersection is a point, so there is an intersection if the point is not an extremity
            // of the edge, and if it is inside the facet
            if constexpr (std::is_same_v<O, Segment3>) {
                if (facet->containsVertex(e->getVertexSrc()) ||
                    facet->containsVertex(e->getVertexDst())) {
                    return false;
                } else {
                    return isInsideWithRayShooting(*ipoint, facet, handle_degree_1_as_ray);
                }
            } else if constexpr (std::is_same_v<O, Ray3>) {
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
            Ray3 r { *r_src->getPoint(), *r_dst->getPoint() };
            return test_o_facet(edge, r);
        }
    }

    Segment3 s { *e_src->getPoint(), *e_dst->getPoint() };
    return test_o_facet(edge, s);
}

// @fixme self-intersection on boundaries does not seem well defined
bool SelfIntersection::hasSelfIntersectingSurface(PolyhedronSPtr polyhedron) {
    CGAL_SS3_ALGO_TRACE("\n> hasSelfIntersectingSurface()");

    if (SelfIntersection::hasSelfIntersectingFacets(polyhedron)) {
        return true;
    }

    // @speed O(nf*ne) algorithm
    for (FacetSPtr facet : polyhedron->facets()) {
        for (EdgeSPtr edge : polyhedron->edges()) {
            if (isEdgeInsideFacet(facet, edge, true /*handle_degree_1_as_ray*/)) {
                CGAL_SS3_ALGO_TRACE("\nPolyhedron has no self-intersecting facets, but the surface is self-intersecting!");
                CGAL_SS3_ALGO_TRACE(facet->toString());
                CGAL_SS3_ALGO_TRACE(edge->toString());
                return true;
            }
        }
    }

    return false;
}

} }
