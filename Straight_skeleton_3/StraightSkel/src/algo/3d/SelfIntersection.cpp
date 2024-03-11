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

Point3SPtr SelfIntersection::intersectEdges(FacetSPtr facet,
        EdgeSPtr edge1, EdgeSPtr edge2, bool handle_deg1_as_ray) {
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
            const double max = std::numeric_limits<double>::max();
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
        if (!reflex) {
            result = KernelWrapper::bisector(
                    KernelWrapper::opposite(plane_in), plane_out);
        } else {
            result = KernelWrapper::bisector(
                    plane_in, KernelWrapper::opposite(plane_out));
        }
    }
    return result;
}

EdgeSPtr SelfIntersection::findNearestEdge(FacetSPtr facet, Point3SPtr point) {
    EdgeSPtr result;
    double dist_min = std::numeric_limits<double>::max();
    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr vertex_src = edge->src(facet);
        VertexSPtr vertex_dst = edge->dst(facet);
        bool point_inside_bounds = true;
        Plane3SPtr bisector_src;
        Plane3SPtr bisector_dst;
        if (vertex_src->degree() > 1) {
            bisector_src = bisector(facet, vertex_src);
        }
        if (vertex_dst->degree() > 1) {
            bisector_dst = bisector(facet, vertex_dst);
        }
        if (bisector_src) {
            if (KernelWrapper::side(bisector_src, point) < 0) {
                point_inside_bounds = false;
            }
        }
        if (bisector_dst) {
            if (KernelWrapper::side(bisector_dst, point) > 0) {
                point_inside_bounds = false;
            }
        }
        if (point_inside_bounds) {
            double dist = KernelWrapper::distance(edge->line(), point);
            if (dist < dist_min) {
                result = edge;
                dist_min = dist;
            }
        }
    }
    return result;
}

bool SelfIntersection::isEdgeInsideFacet(FacetSPtr facet, EdgeSPtr edge, bool handle_deg1_as_ray) {
    bool result = false;
    if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
        return false;
    }
    if (facet->containsVertex(edge->getVertexSrc()) ||
            facet->containsVertex(edge->getVertexDst())) {
        return false;
    }
    Line3SPtr line = edge->line();
    Point3SPtr point = KernelWrapper::intersection(facet->plane(), line);
    if (point) {
        Point3SPtr p_src = edge->getVertexSrc()->getPoint();
        Point3SPtr p_dst = edge->getVertexDst()->getPoint();
        if (handle_deg1_as_ray) {
            // 0 * inf = nan
            // 0 * max = 0
            const double max = std::numeric_limits<double>::max();
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
