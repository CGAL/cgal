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
#include <limits>
#include <list>

namespace algo { namespace _3d {

// @todo all of this should be a predicate
Point3SPtr SelfIntersection::intersectEdges(FacetSPtr facet,
        EdgeSPtr edge1, EdgeSPtr edge2, bool handle_deg1_as_ray)
{
    Point3SPtr result = Point3SPtr();

    // @fixme? full overlaps?
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
        Plane3SPtr plane1 = KernelFactory::createPlane3(p_11, p_12, p_13);
        Point3SPtr p_21 = edge2->getVertexSrc()->getPoint();
        Point3SPtr p_22 = edge2->getVertexDst()->getPoint();
        Point3SPtr p_23 = KernelFactory::createPoint3((*p_21) + (*normal));
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

        if (handle_deg1_as_ray) {
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
    return result;
}

bool SelfIntersection::isSelfIntersectingFacet(FacetSPtr facet) {
    bool result = false;
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
                result = true;
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
            result++;
        }
    }
    return result;
}

// #define CGAL_SS3_OLD_CODE_FIND_NEAREST_EDGE_CODE
#ifdef CGAL_SS3_OLD_CODE_FIND_NEAREST_EDGE_CODE
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
#else // CGAL_SS3_OLD_CODE_FIND_NEAREST_EDGE_CODE
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
# error "This function is not compatible with the old kernel
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

        // @todo only compute this when necessary and it should be a predicate
        CGAL::FT sq_dist_in = KernelWrapper::squared_distance(edge_in->line(), point);
        CGAL::FT sq_dist_out = KernelWrapper::squared_distance(edge_out->line(), point);
        std::cout << "sq_dist_in = " << sq_dist_in << std::endl;
        std::cout << "sq_dist_out = " << sq_dist_out << std::endl;

        if(or_in == CGAL::POSITIVE) {
            if(or_out == CGAL::POSITIVE) {
                // upper quadrant
                result = CGAL::compare(sq_dist_out, sq_dist_in);
            } else { // or_out != CGAL::POSITIVE
                // in quadrant (left in convex, right in reflex)
                if(*point == *p_mid) {
                  result = CGAL::ON_ORIENTED_BOUNDARY;
                } else {
                  result = CGAL::ON_NEGATIVE_SIDE;
                }
            }
        } else { // or_in != CGAL::POSITIVE
            if(or_out == CGAL::POSITIVE) {
                // out quadrant (left in reflex, right in convex)
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
    } else {
        CGAL_unreachable();
    }

    return result;
}
#endif // CGAL_SS3_OLD_CODE_FIND_NEAREST_EDGE_CODE

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
#ifdef CGAL_SS3_OLD_CODE_FIND_NEAREST_EDGE_CODE
            Plane3SPtr bisector_src = bisector(facet, vertex_src);
            if (bisector_src) {
                if (KernelWrapper::side(bisector_src, point) < 0) {
                    point_inside_bounds = false;
                }
            }

            // if(!point_inside_bounds) {
            //    CGAL_assertion(SideOfBisector(facet, vertex_src, point) == CGAL::ON_NEGATIVE_SIDE);
            // } else {
            //    CGAL_assertion(SideOfBisector(facet, vertex_src, point) != CGAL::ON_NEGATIVE_SIDE);
            // }
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

                // if(!point_inside_bounds) {
                //    CGAL_assertion(SideOfBisector(facet, vertex_src, point) == CGAL::ON_POSITIVE_SIDE);
                // } else {
                //    CGAL_assertion(SideOfBisector(facet, vertex_src, point) != CGAL::ON_POSITIVE_SIDE);
                // }
#else
                CGAL::Sign res = SideOfBisector(facet, vertex_dst, point);
                if (res == CGAL::ON_POSITIVE_SIDE) {
                    point_inside_bounds = false;
                }
#endif
            }
        }

        if (point_inside_bounds) {
            CGAL::FT sq_dist = KernelWrapper::squared_distance(edge->line(), point);
            if (sq_dist < sq_dist_min) { // @fixme construction
                result = edge;
                sq_dist_min = sq_dist;
            }
        }
    }
    return result;
}

bool SelfIntersection::isEdgeInsideFacet(FacetSPtr facet, EdgeSPtr edge, bool handle_deg1_as_ray)
{
    bool result = false;
    if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
        return false;
    }

    // @fixme? isn't it possible to have an edge coplanar with the facet?
    if (facet->containsVertex(edge->getVertexSrc()) ||
            facet->containsVertex(edge->getVertexDst())) {
        return false;
    }

    Line3SPtr line = edge->line();

    // @fixme? what if the intersection isn't a point?
    Point3SPtr point = KernelWrapper::intersection(facet->plane(), line);
    if (point) {
        Point3SPtr p_src = edge->getVertexSrc()->getPoint();
        Point3SPtr p_dst = edge->getVertexDst()->getPoint();
        if (handle_deg1_as_ray) {
            // 0 * inf = nan
            // 0 * max = 0
            const CGAL::FT max = std::numeric_limits<double>::max();
            if (edge->getVertexSrc()->degree() == 1) {
                p_src = KernelFactory::createPoint3(*p_src +
                        (*p_src-*p_dst) * max);
            }
            if (edge->getVertexDst()->degree() == 1) {
                p_dst = KernelFactory::createPoint3(*p_dst +
                        (*p_dst-*p_src) * max);
            }
        }
        if (KernelWrapper::isInside(point, p_src, p_dst)) {
            Vector3SPtr normal = KernelFactory::createVector3(facet->plane());
            Line3SPtr line_point = KernelFactory::createLine3(point, normal);
            EdgeSPtr edge_nearest = findNearestEdge(facet, point);
            Line3SPtr line_nearest = edge_nearest->line();
            if (edge_nearest->getFacetL() == facet) {
                // @fixme the point could be aligned with the nearest edge
                // (but currently it's fine because the closest edge is found using
                // bisectors)
                if (KernelWrapper::orientation(line_nearest, line_point) >= 0) {
                    result = true;
                }
            } else if (edge_nearest->getFacetR() == facet) {
                if (KernelWrapper::orientation(line_nearest, line_point) <= 0) {
                    result = true;
                }
            }
        }
    }
    return result;
}

bool SelfIntersection::hasSelfIntersectingSurface(PolyhedronSPtr polyhedron) {
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
                if (isEdgeInsideFacet(facet, edge, true)) {
                    DEBUG_PRINT("Polyhedron has no self-intersecting facets, but the surface is self-intersecting.");
                    result = true;
                    break;
                }
            }
            if (result) {
                break;
            }
        }
    }
    return result;
}

} }
