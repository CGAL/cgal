/**
 * @file   algo/3d/PolyhedronIntersection.cpp
 * @author Gernot Walzl
 * @date   2012-09-26
 */

#include "algo/3d/PolyhedronIntersection.h"

#include "algo/3d/KernelWrapper.h"
#include "data/2d/Edge.h"
#include "data/2d/Polygon.h"
#include "data/2d/skel/SkelEdgeData.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/skel/SkelEdgeData.h"
#include <cmath>
#include <limits>
#include <list>

namespace algo { namespace _3d {

PolyhedronIntersection::PolyhedronIntersection() {
    // intentionally does nothing.
}

PolyhedronIntersection::~PolyhedronIntersection() {
    // intentionally does nothing.
}

Point3SPtr PolyhedronIntersection::intersect(EdgeSPtr edge, Plane3SPtr plane) {
    Point3SPtr result = Point3SPtr();
    Point3SPtr p = KernelWrapper::intersection(plane, edge->line());
    if (p) {
        if (KernelWrapper::isInside(p,
                edge->getVertexSrc()->getPoint(),
                edge->getVertexDst()->getPoint())) {
            result = p;
        }
    }
    return result;
}

EdgeSPtr PolyhedronIntersection::findDst(FacetSPtr facet, EdgeSPtr edge_src,
        Plane3SPtr plane) {
    EdgeSPtr result = EdgeSPtr();

    Point3SPtr p_src = intersect(edge_src, plane);
    if (!p_src) {
        return result;
    }

    Vector3SPtr normal_plane = KernelFactory::createVector3(plane);
    Vector3SPtr normal_poly = KernelFactory::createVector3(facet->plane());
    Vector3SPtr dir_req = KernelWrapper::cross(normal_plane, normal_poly);

    double dist_min = std::numeric_limits<double>::max();
    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (edge == edge_src) {
            continue;
        }
        Point3SPtr p_dst = intersect(edge, plane);
        if (p_dst) {
            Vector3SPtr dir = KernelFactory::createVector3((*p_dst)-(*p_src));
            //if (KernelWrapper::angle(dir_req, dir) < M_PI/2.0) {
            if ((*dir_req * *dir) > 0.0) {
                double dist = KernelWrapper::distance(p_src, p_dst);
                if (dist < dist_min) {
                    result = edge;
                    dist_min = dist;
                }
            }
        }
    }
    return result;
}

FacetSPtr PolyhedronIntersection::intersect(PolyhedronSPtr polyhedron,
        Plane3SPtr plane) {
    FacetSPtr result = Facet::create();
    result->setPlane(plane);

    EdgeSPtr edge_begin;
    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        Point3SPtr p = intersect(edge, plane);
        if (p) {
            edge_begin = edge;
            break;
        }
    }

    FacetSPtr facet;
    if (KernelWrapper::side(plane, edge_begin->getVertexSrc()->getPoint()) > 0) {
        facet = edge_begin->getFacetL();
    } else {
        facet = edge_begin->getFacetR();
    }

    VertexSPtr vertex_begin;
    VertexSPtr vertex_src;
    VertexSPtr vertex_dst;
    EdgeSPtr edge_src;
    EdgeSPtr edge_dst;
    while (edge_src != edge_begin) {
        if (!edge_src) {
            edge_src = edge_begin;
            Point3SPtr p_src = intersect(edge_src, plane);
            vertex_src = Vertex::create(p_src);
            result->addVertex(vertex_src);
            vertex_begin = vertex_src;
        }
        edge_dst = findDst(facet, edge_src, plane);
        if (edge_dst == edge_begin) {
            vertex_dst = vertex_begin;
        } else {
            Point3SPtr p_dst = intersect(edge_dst, plane);
            vertex_dst = Vertex::create(p_dst);
            result->addVertex(vertex_dst);
        }

        EdgeSPtr edge = Edge::create(vertex_src, vertex_dst);
        edge->setFacetL(result);
        skel::SkelEdgeDataSPtr data = skel::SkelEdgeData::create(edge);
        data->setFacetOrigin(facet);
        result->addEdge(edge);

        vertex_src = vertex_dst;
        edge_src = edge_dst;
        facet = edge_src->other(facet);
    }
    return result;
}

data::_2d::PolygonSPtr PolyhedronIntersection::toWeighted2d(FacetSPtr facet) {
    data::_2d::PolygonSPtr result = facet->toPolygon();
    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    std::list<data::_2d::EdgeSPtr>::iterator it_e2 = result->edges().begin();
    while (it_e != facet->edges().end() && it_e2 != result->edges().end()) {
        EdgeSPtr edge = *it_e++;
        data::_2d::EdgeSPtr edge2 = *it_e2++;
        skel::SkelEdgeDataSPtr data = std::dynamic_pointer_cast<skel::SkelEdgeData>(
                edge->getData());
        double angle = KernelWrapper::angle(
                facet->plane(), data->getFacetOrigin()->plane());
        double speed = 1.0/sin(angle);
        data::_2d::skel::SkelEdgeDataSPtr data2 =
                data::_2d::skel::SkelEdgeData::create(edge2);
        data2->setSpeed(speed);
    }
    return result;
}

} }

