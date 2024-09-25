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

#include <limits>
#include <list>

#define CGAL_SS3_EXIT_ASAP

namespace algo { namespace _3d {

// whether the edges intersect in their interior, with edges possibly being rays
// @todo this must be rewritten to be predicates
bool SelfIntersection::intersectEdges(FacetSPtr facet,
                                      EdgeSPtr edge1,
                                      EdgeSPtr edge2,
                                      bool handle_degree_1_as_ray)
{
    CGAL_assertion(edge1 != edge2);
    CGAL_assertion(edge1->getVertexSrc()->degree() != 1 || edge1->getVertexDst()->degree() != 1);
    CGAL_assertion(edge2->getVertexSrc()->degree() != 1 || edge2->getVertexDst()->degree() != 1);

// #define CGAL_SS3_OLD_CODE_INTERSECT_EDGES
#ifdef CGAL_SS3_OLD_CODE_INTERSECT_EDGES
    Point3SPtr result = Point3SPtr();

    if (edge1->getVertexSrc() == edge2->getVertexSrc() ||
            edge1->getVertexSrc() == edge2->getVertexDst()) {
        result = edge1->getVertexSrc()->getPoint();
    } else if (edge1->getVertexDst() == edge2->getVertexSrc() ||
            edge1->getVertexDst() == edge2->getVertexDst()) {
        result = edge1->getVertexDst()->getPoint();
    } else {
        Vector3SPtr normal = KernelFactory::createVector3(facet->plane());
        Point3SPtr p_11 = edge1->getVertexSrc()->getPoint();
        Point3SPtr p_12 = edge1->getVertexDst()->getPoint();
        Point3SPtr p_13 = KernelFactory::createPoint3((*p_11) + (*normal));
        CGAL_assertion(!CGAL::collinear(*p_11, *p_12, *p_13));
        Plane3SPtr plane1 = KernelFactory::createPlane3(p_11, p_12, p_13);
        Point3SPtr p_21 = edge2->getVertexSrc()->getPoint();
        Point3SPtr p_22 = edge2->getVertexDst()->getPoint();
        Point3SPtr p_23 = KernelFactory::createPoint3((*p_21) + (*normal));
        CGAL_assertion(!CGAL::collinear(*p_21, *p_22, *p_23));
        Plane3SPtr plane2 = KernelFactory::createPlane3(p_21, p_22, p_23);
        Point3SPtr p_intersection = KernelWrapper::intersection(plane1, plane2,
                facet->plane());

        // fix rounding errors (they cause problems in isInside)
        for (unsigned int i = 0; i < 3; i++) {
            if ((*p_11)[i] == (*p_12)[i]) {
                p_intersection =
                        KernelWrapper::replaceCoord(p_intersection, p_11, i);
            } else if ((*p_21)[i] == (*p_22)[i]) {
                p_intersection =
                        KernelWrapper::replaceCoord(p_intersection, p_21, i);
            }
        }

        if (handle_degree_1_as_ray) {
            // 0 * inf = nan
            // const double infinity = std::numeric_limits<double>::infinity();
            // 0 * max = 0
            const CGAL::FT max = std::numeric_limits<double>::max(); // do not put FT
            if (edge1->getVertexSrc()->degree() == 1) {
                p_11 = KernelFactory::createPoint3(*p_11 +
                        (*p_11-*p_12) * max);
            }
            if (edge1->getVertexDst()->degree() == 1) {
                p_12 = KernelFactory::createPoint3(*p_12 +
                        (*p_12-*p_11) * max);
            }
            if (edge2->getVertexSrc()->degree() == 1) {
                p_21 = KernelFactory::createPoint3(*p_21 +
                        (*p_21-*p_22) * max);
            }
            if (edge2->getVertexDst()->degree() == 1) {
                p_22 = KernelFactory::createPoint3(*p_22 +
                        (*p_22-*p_21) * max);
            }
        }

        if (KernelWrapper::isInside(p_intersection, p_11, p_12) &&
                KernelWrapper::isInside(p_intersection, p_21, p_22)) {
            result = p_intersection;
        }
    }

    return (result != Point3SPtr());
#else
    // if there are intersections, tolerate it only if it's an extremity
    auto treat_seg_seg = [](const Segment3& s1, const Segment3& s2) -> bool
    {
        CGAL::Object obj = CGAL::intersection(s1, s2);
        if (!obj) {
          return false;
        } else if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
            return (*ipoint != s1.source() && *ipoint != s1.target() &&
                    *ipoint != s2.source() && *ipoint != s2.target());
        } else {
            return true; // intersection is a segment
        }
    };

    auto treat_seg_ray = [](const Segment3& s, const Ray3& r) -> bool
    {
        CGAL::Object obj = CGAL::intersection(s, r);
        if (!obj) {
          return false;
        } else if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
            return (*ipoint != s.source() && *ipoint != s.target() && *ipoint != r.source());
        } else {
            return true; // intersection is a segment
        }
    };

    auto treat_ray_ray = [](const Ray3& r1, const Ray3& r2) -> bool
    {
        CGAL::Object obj = CGAL::intersection(r1, r2);
        if (!obj) {
          return false;
        } else if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
            return (*ipoint != r1.source() && *ipoint != r2.source());
        } else {
            return true; // intersection is a segment
        }
    };

    Point3SPtr p_11 = edge1->getVertexSrc()->getPoint();
    Point3SPtr p_12 = edge1->getVertexDst()->getPoint();
    Point3SPtr p_21 = edge2->getVertexSrc()->getPoint();
    Point3SPtr p_22 = edge2->getVertexDst()->getPoint();

    // skip if one is degenerate: if there is a real intersection,
    // we'll meet it through a non-degenerate combination
    if(*p_11 == *p_12 || *p_21 == *p_22) {
        return false;
    }

    if (handle_degree_1_as_ray) {
        if (edge1->getVertexSrc()->degree() == 1) {
            if (edge2->getVertexSrc()->degree() == 1) {
                return treat_ray_ray(Ray3{*p_12,*p_11}, Ray3{*p_22,*p_21});
            } else if(edge2->getVertexDst()->degree() == 1) {
                return treat_ray_ray(Ray3{*p_12,*p_11}, Ray3{*p_21,*p_22});
            } else {
                return treat_seg_ray(Segment3{*p_21,*p_22}, Ray3{*p_12,*p_11});
            }
        } else if (edge1->getVertexDst()->degree() == 1) {
            if (edge2->getVertexSrc()->degree() == 1) {
                return treat_ray_ray(Ray3{*p_11,*p_12}, Ray3{*p_22,*p_21});
            } else if(edge2->getVertexDst()->degree() == 1) {
                return treat_ray_ray(Ray3{*p_11,*p_12}, Ray3{*p_21,*p_22});
            } else {
                return treat_seg_ray(Segment3{*p_21,*p_22}, Ray3{*p_11,*p_12});
            }
        } else { // no degree 1 in edge1
            if (edge2->getVertexSrc()->degree() == 1) {
                return treat_seg_ray(Segment3{*p_11,*p_12}, Ray3{*p_22,*p_21});
            } else if(edge2->getVertexDst()->degree() == 1) {
                return treat_seg_ray(Segment3{*p_11,*p_12}, Ray3{*p_21,*p_22});
            }
            // neither edge having degree 1 vertices is below
        }
    }

    return treat_seg_seg(Segment3{*p_11, *p_12}, Segment3{*p_21,*p_22});
#endif
}

// @todo currently has square complexity, it could be a sweep like in "is_simple_polygon_2"
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
            if (edge1->getVertexSrc() == edge2->getVertexSrc() ||
                    edge1->getVertexSrc() == edge2->getVertexDst() ||
                    edge1->getVertexDst() == edge2->getVertexSrc() ||
                    edge1->getVertexDst() == edge2->getVertexDst()) {
                continue;
            }
            if (intersectEdges(facet, edge1, edge2, true)) {
                std::cout << "edges intersect within the face:" << std::endl;
                std::cout << facet->toString() << std::endl;
                std::cout << edge1->toString() << std::endl;
                std::cout << edge2->toString() << std::endl;
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
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
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

Plane3SPtr SelfIntersection::bisector(FacetSPtr facet, VertexSPtr vertex) {
    Plane3SPtr result;
    EdgeSPtr edge_in;
    EdgeSPtr edge_out;
    std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
    while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
            if (edge->src(facet) == vertex) {
                edge_out = edge;
            } else if (edge->dst(facet) == vertex) {
                edge_in = edge;
            }
        }
    }
    if (edge_in && edge_out) {
        Point3SPtr point = vertex->getPoint();
        Point3SPtr p_in = edge_in->src(facet)->getPoint();
        Point3SPtr p_out = edge_out->dst(facet)->getPoint();
        Vector3SPtr normal = KernelFactory::createVector3(facet->plane());
        Point3SPtr p_normal = KernelFactory::createPoint3(*point + *normal);
        Plane3SPtr plane_in = KernelFactory::createPlane3(p_in, point, p_normal);
        Plane3SPtr plane_out = KernelFactory::createPlane3(point, p_out, p_normal);
        bool reflex = false;
        if (KernelWrapper::side(plane_in, p_out) > 0) {
            reflex = true;
        }

        // p_out is on the positive side of the bisector
        if (!reflex) {
            result = KernelWrapper::bisector(
                    KernelWrapper::opposite(plane_in), plane_out);
        } else {
            result = KernelWrapper::bisector(
                    plane_in, KernelWrapper::opposite(plane_out));
        }

        CGAL_postcondition(KernelWrapper::side(result, p_out) >= 0);
    }
    return result;
}

// Returns:
// - ON_POSITIVE_SIDE  if the query is on the positive side (see below for definition)
//   of the planar bisector of the two edges incident to 'vertex'.
// - ON_ORIENTED_BOUNDARY if the query is on the bisector
// - ON_NEGATIVE_SIDE if the query is on the negative side of the bisector
//
// The positive side of the bisector is the side that has p_out
// in the following (facet seen from above):
//
//                         | <-- planar bisector (positive side is here)
//                         |
//           edge_in       |          edge_out
//  p_in --------------> vertex -------------------> p_out
//                         |
//                         |
//
// \pre `point` is coplanar with the edges incident to `vertex`
CGAL::Sign SelfIntersection::SideOfBisector(FacetSPtr facet, VertexSPtr vertex, Point3SPtr point) {
#ifndef USE_CGAL
# error "This function is not compatible with the old kernel"
#endif

    CGAL_precondition(vertex->degree() > 1);

    CGAL::Sign result;
    EdgeSPtr edge_in;
    EdgeSPtr edge_out;

    std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
    while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
            if (edge->src(facet) == vertex) {
                edge_out = edge;
            } else if (edge->dst(facet) == vertex) {
                edge_in = edge;
            }
        }
    }

    if (edge_in && edge_out) {
        // @fixme plenty of constructions below
        Point3SPtr p_mid = vertex->getPoint();
        Point3SPtr p_in = edge_in->src(facet)->getPoint();
        Point3SPtr p_out = edge_out->dst(facet)->getPoint();
        Vector3SPtr normal = KernelFactory::createVector3(facet->plane());
        Point3SPtr p_normal = KernelFactory::createPoint3(*point + *normal);

        CGAL::Orientation or_in = CGAL::orientation(*p_mid, *p_normal, *p_in, *point);
        CGAL::Orientation or_out = CGAL::orientation(*p_mid, *p_out, *p_normal, *point);

        // positive = right turn (reflex), null = collinear, negative = left turn (convex)
        CGAL::Orientation or_edge = CGAL::orientation(*p_mid, *p_normal, *p_in, *p_out);

        // @todo only compute this when necessary and it should be a predicate
        // Note that there are no weights involved here: we are checking for self-intersections
        // in a fixed polyhedron
        CGAL::FT sq_dist_in = KernelWrapper::squared_distance(edge_in->line(), point);
        CGAL::FT sq_dist_out = KernelWrapper::squared_distance(edge_out->line(), point);

        // std::cout << "p_in = " << *p_in << std::endl;
        // std::cout << "p_mid = " << *p_mid << std::endl;
        // std::cout << "p_out = " << *p_out << std::endl;
        // std::cout << "query = " << *point << std::endl;
        // std::cout << "or_in = " << or_in << std::endl;
        // std::cout << "or_out = " << or_out << std::endl;
        // std::cout << "or_edge = " << or_edge << std::endl;
        // std::cout << "sq_dist_in = " << sq_dist_in << std::endl;
        // std::cout << "sq_dist_out = " << sq_dist_out << std::endl;

        // that's for the left turn (convex); for reflex, it's the opposite
        if(or_in == CGAL::POSITIVE) {
            if(or_out == CGAL::POSITIVE) {
                // upper quadrant
                result = CGAL::compare(sq_dist_out, sq_dist_in);
            } else { // or_out != CGAL::POSITIVE
                // right quadrant
                if(*point == *p_mid) {
                  result = CGAL::ON_ORIENTED_BOUNDARY;
                } else {
                  result = CGAL::ON_NEGATIVE_SIDE;
                }
            }
        } else { // or_in != CGAL::POSITIVE
            if(or_out == CGAL::POSITIVE) {
                // left quadrant
                if(*point == *p_mid) {
                    result = CGAL::ON_ORIENTED_BOUNDARY;
                } else {
                    result = CGAL::ON_POSITIVE_SIDE;
                }
            } else { // or_out != CGAL::POSITIVE
                // bottom quadrant
                result = CGAL::compare(sq_dist_in, sq_dist_out);
            }
        }

        // opposite for reflex
        if (or_edge == CGAL::POSITIVE /*reflex*/) {
          result = (result == CGAL::POSITIVE) ? CGAL::NEGATIVE : CGAL::POSITIVE;
        }

    } else {
        CGAL_unreachable();
    }

    return result;
}

// @todo could be improved by storing the previous (query, edge) position
// --> for consecutive edges only!
EdgeSPtr SelfIntersection::findNearestEdge(FacetSPtr facet, Point3SPtr point) {
    EdgeSPtr result;
    CGAL::FT sq_dist_min = std::numeric_limits<double>::max(); // do not put FT
    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr vertex_src = edge->src(facet);
        VertexSPtr vertex_dst = edge->dst(facet);

        bool point_inside_bounds = true;

        if (vertex_src->degree() > 1) {
#define CGAL_SS3_OLD_CODE_FIND_NEAREST_EDGE_CODE
#ifdef CGAL_SS3_OLD_CODE_FIND_NEAREST_EDGE_CODE
            Plane3SPtr bisector_src = bisector(facet, vertex_src);
            if (bisector_src) {
                if (KernelWrapper::side(bisector_src, point) < 0) {
                    point_inside_bounds = false;
                }
            }

            CGAL_assertion_code(
                auto sob_res = SideOfBisector(facet, vertex_src, point);
                std::cout << "SideOfBisector(facet, vertex_src, point) = " << sob_res << std::endl;
                std::cout << "point_inside_bounds = " << point_inside_bounds << std::endl;
            )
            if(point_inside_bounds) {
               CGAL_assertion(sob_res != CGAL::ON_NEGATIVE_SIDE);
            } else {
               CGAL_assertion(sob_res == CGAL::ON_NEGATIVE_SIDE);
            }
#else
            CGAL::Sign res = SideOfBisector(facet, vertex_src, point);
            if (res == CGAL::ON_NEGATIVE_SIDE) {
                point_inside_bounds = false;
            }
#endif
        }

        if(point_inside_bounds) {
            if (vertex_dst->degree() > 1) {
#ifdef CGAL_SS3_OLD_CODE_FIND_NEAREST_EDGE_CODE
                Plane3SPtr bisector_dst = bisector(facet, vertex_dst);
                if (bisector_dst) {
                    if (KernelWrapper::side(bisector_dst, point) > 0) {
                        point_inside_bounds = false;
                    }
                }

                CGAL_assertion_code(
                    auto sob_res = SideOfBisector(facet, vertex_dst, point);
                    std::cout << "point_inside_bounds = " << point_inside_bounds << std::endl;
                    std::cout << "SideOfBisector(facet, vertex_dst, point) = " << sob_res << std::endl;
                )
                if(point_inside_bounds) {
                   CGAL_assertion(sob_res != CGAL::ON_POSITIVE_SIDE);
                } else {
                   CGAL_assertion(sob_res == CGAL::ON_POSITIVE_SIDE);
                }
#else
                CGAL::Sign res = SideOfBisector(facet, vertex_dst, point);
                if (res == CGAL::ON_POSITIVE_SIDE) {
                    point_inside_bounds = false;
                }
#endif
            }
        }

        if (point_inside_bounds) {
             // @fixme construction
            CGAL::FT sq_dist = KernelWrapper::squared_distance(edge->line(), point);
            if (sq_dist < sq_dist_min) {
                result = edge;
                sq_dist_min = sq_dist;
            }
        }
    }
    return result;
}

// @todo construction galore...
bool SelfIntersection::isInsideWithRayShooting(Point3SPtr point,
                                               FacetSPtr facet)
{
    std::cout << "\n> isInsideWithRayShooting()" << std::endl;
    std::cout << facet->toString() << std::endl;

    Plane3SPtr pl = facet->plane();
    Vector3SPtr normal = KernelFactory::createVector3(pl);

    // shoot random rays till something is hit
    // essential to this: we know the facet does not self-intersect
    CGAL::Random rng(0);
    CGAL::Random_points_on_sphere_3<Point3> random_point_on_sphere(1, rng);

    // some initial rays that we know are in the plane, and know will hit something
    std::vector<Point3> candidate_ray_targets;

    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        Point3SPtr p_src = edge->getVertexSrc()->getPoint();
        Point3SPtr p_dst = edge->getVertexDst()->getPoint();

        // std::cout << " p_src = " << *p_src << std::endl;
        // std::cout << " p_dst = " << *p_dst << std::endl;

        // @todo this only stands for EPECK and faces merged with zero tolerance
        CGAL_assertion(pl->has_on(*point));
        CGAL_assertion(pl->has_on(*p_src));
        CGAL_assertion(pl->has_on(*p_dst));
        CGAL_assertion(CGAL::scalar_product(Vector3(*p_src, *p_dst), *normal) == 0);

        candidate_ray_targets.push_back(CGAL::midpoint(*p_src, *p_dst));
    }

    // @fixme handle the case of an half plane face with the query point on the boundary

    for(;;) {
        Ray3 shooting_ray;
        if (!candidate_ray_targets.empty()) {
            shooting_ray = Ray3(*point, candidate_ray_targets.back());
            candidate_ray_targets.pop_back();
        } else {
            Point3 rnd_p = *random_point_on_sphere++;
            Point3 target_p = *point + Vector3(CGAL::ORIGIN, rnd_p);
            Point3 proj_p = pl->projection(target_p);
            if(proj_p == *point) {
                continue;
            }
            shooting_ray = Ray3(*point, proj_p);
        }

        // std::cout << "shooting_ray = " << shooting_ray.point(0) << " " << shooting_ray.point(1) << std::endl;

        CGAL_assertion(shooting_ray.point(0) == *point);
        CGAL_assertion(pl->has_on(shooting_ray.point(0)));
        CGAL_assertion(pl->has_on(shooting_ray.point(1)));

        CGAL::FT sq_dist_to_closest = std::numeric_limits<double>::max();
        EdgeSPtr closest_edge;

        auto treat_edge = [&](EdgeSPtr target_edge, const auto& edge_geometry) -> void
        {
            CGAL::Object obj = CGAL::intersection(shooting_ray, edge_geometry);
            if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
                CGAL::FT sqd = CGAL::squared_distance(*point, *ipoint);
                // std::cout << "intersects & sq_dst: " << sqd << std::endl;
                CGAL::Comparison_result res = CGAL::compare(sqd, sq_dist_to_closest);
                if (res == CGAL::SMALLER) {
                    sq_dist_to_closest = sqd;
                    closest_edge = target_edge;
                }
            }
        };

        // @todo brute force finding the closest, but faces are small so...
        std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            // std::cout << "consider edge: " << edge->toString() << std::endl;
            VertexSPtr v_src = edge->src(facet);
            VertexSPtr v_dst = edge->dst(facet);
            Point3SPtr p_src = v_src->getPoint();
            Point3SPtr p_dst = v_dst->getPoint();

            // We'll use a (coplanar) orientation check with the edge to determine
            // whether we are on the 'outside' or 'inside'
            if (CGAL::collinear(*p_src, *p_dst, *point)) {
                if (CGAL::collinear_are_ordered_along_line(*p_src, *point, *p_dst)) {
                    return true; // on edge
                } else {
                    std::cout << " skipping because collinear" << std::endl;
                }
                continue;
            }

            if (v_src->degree() == 1) {
                if (v_dst->degree() == 1) {
                    // std::cout << "L3" << std::endl;
                    treat_edge(edge, Line3(*p_src, *p_dst));
                } else {
                    // std::cout << "R3 from DST " << *p_dst << std::endl;
                    treat_edge(edge, Ray3(*p_dst, *p_src));
                }
            } else if (v_dst->degree() == 1) {
                // std::cout << "R3 from SRC " << *p_src << std::endl;
                treat_edge(edge, Ray3(*p_src, *p_dst));
            } else {
                // std::cout << "S3" << std::endl;
                treat_edge(edge, Segment3(*p_src, *p_dst));
            }
        }

        // if (closest_edge)
        //     std::cout << "closest_edge = " << closest_edge->toString() << std::endl;
        // else
        //     std::cout << "closest_edge = NONE" << std::endl;

        // Being in the cone where the distance is the same to multiple edges
        // makes the orientation test not usable
        if (closest_edge == EdgeSPtr()) {
            continue; // try another ray; we'll hit something eventually!
        }

        Point3SPtr p_src = closest_edge->src(facet)->getPoint();
        Point3SPtr p_dst = closest_edge->dst(facet)->getPoint();

        // std::cout << "p_src = " << *p_src << std::endl;
        // std::cout << "p_src + normal = " << *p_src + *normal << std::endl;
        // std::cout << "p_dst = " << *p_dst << std::endl;
        // std::cout << "point = " << *point << std::endl;
        CGAL_assertion(!CGAL::collinear(*p_src, *p_src + *normal, *p_dst));
        CGAL_assertion(CGAL::scalar_product(Vector3(*p_src, *p_src + *normal), Vector3(*p_src, *p_dst)) == 0);

        CGAL::Orientation o = CGAL::orientation(*p_src, *p_src + *normal, *p_dst, *point);
        // std::cout << "Orientation = " << o << std::endl;

        return (o != CGAL::NEGATIVE);
    }

    CGAL_unreachable();
    return false;
}

bool SelfIntersection::isEdgeInsideFacet(FacetSPtr facet,
                                         EdgeSPtr edge,
                                         bool handle_degree_1_as_ray)
{
    // std::cout << "\n> isEdgeInsideFacet()" << std::endl;

    bool result = false;

    if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
        return false;
    }

    // @fixme? isn't it possible to have an edge coplanar with the facet?
    if (facet->containsVertex(edge->getVertexSrc()) ||
            facet->containsVertex(edge->getVertexDst())) {
        return false;
    }

    if (edge->getVertexSrc()->getPoint() == edge->getVertexDst()->getPoint()) {
        return false;
    }

    // @todo could do with a segment or ray (if degree 1) direclty
    Line3SPtr line = edge->line();

    Point3SPtr point = KernelWrapper::intersection(facet->plane(), line);
    CGAL_assertion(bool(point)); // @fixme? what if the intersection isn't a point?
    // std::cout << "Intersection point " << *point << std::endl;

    if (point) {
        Point3SPtr p_src = edge->getVertexSrc()->getPoint();
        Point3SPtr p_dst = edge->getVertexDst()->getPoint();
        if (handle_degree_1_as_ray) {
            // 0 * inf = nan
            // 0 * max = 0
            const CGAL::FT max = std::numeric_limits<double>::max();
            if (edge->getVertexSrc()->degree() == 1) {
                p_src = KernelFactory::createPoint3(*p_src + (*p_src-*p_dst) * max);
            }
            if (edge->getVertexDst()->degree() == 1) {
                p_dst = KernelFactory::createPoint3(*p_dst + (*p_dst-*p_src) * max);
            }
        }

        if (KernelWrapper::isInside(point, p_src, p_dst)) {

// #define CGAL_SS3_USE_OLD_CODE_FOR_INSIDE_OUT_CHECKS
#ifdef CGAL_SS3_USE_OLD_CODE_FOR_INSIDE_OUT_CHECKS
            // This old code is broken because the findNearestEdge might not return the correct
            // edge.

            // !! WARNING!! if you use this stuff again, you can't call the intersection function
            // with degenerate edges (so, no checking for self intersections after event handling
            // for example)

            Vector3SPtr normal = KernelFactory::createVector3(facet->plane());
            Line3SPtr line_point = KernelFactory::createLine3(point, normal);
            EdgeSPtr edge_nearest = findNearestEdge(facet, point);
            Line3SPtr line_nearest = edge_nearest->line();
            std::cout << "Nearest edge " << edge_nearest->toString() << std::endl;

            if (edge_nearest->getFacetL() == facet) {
                if (KernelWrapper::orientation(line_nearest, line_point) >= 0) {
# ifdef CGAL_SS3_EXIT_ASAP
                    return true;
# else
                    result = true;
# endif
                }
            } else if (edge_nearest->getFacetR() == facet) {
                if (KernelWrapper::orientation(line_nearest, line_point) <= 0) {
# ifdef CGAL_SS3_EXIT_ASAP
                    return true;
# else
                    result = true;
# endif
                }
            }
#else
            if (isInsideWithRayShooting(point, facet)) {
//                CGAL_assertion(result); // @tmp
# ifdef CGAL_SS3_EXIT_ASAP
                return true;
# else
                result = true;
# endif
            } else {
//              CGAL_assertion(!result); // @tmp
            }
#endif
        }
    }

    return result;
}

// @fixme self-intersection on boundaries doesn't seem solidly defined
bool SelfIntersection::hasSelfIntersectingSurface(PolyhedronSPtr polyhedron) {
    // std::cout << "\n> hasSelfIntersectingSurface()" << std::endl;

    bool result = false;
    if (SelfIntersection::hasSelfIntersectingFacets(polyhedron)) {
        result = true;
    } else {
        std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
        while (it_f != polyhedron->facets().end()) {
            FacetSPtr facet = *it_f++;
            std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
            while (it_e != polyhedron->edges().end()) {
                EdgeSPtr edge = *it_e++;
                if (isEdgeInsideFacet(facet, edge, true /*handle_degree_1_as_ray*/)) {
                    std::cout << "\nPolyhedron has no self-intersecting facets, but the surface is self-intersecting!" << std::endl;
                    std::cout << facet->toString() << std::endl;
                    std::cout << edge->toString() << std::endl;
                    result = true;
                    break;
                }
            }

            if (result) { // @fixme just return early...
                break;
            }
        }
    }

    return result;
}

} }
