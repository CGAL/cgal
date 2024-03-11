/**
 * @file   algo/3d/WeightVertexSplitter.cpp
 * @author Gernot Walzl
 * @date   2012-10-03
 */

#include "algo/3d/WeightVertexSplitter.h"

#include "algo/Controller.h"
#include "algo/2d/ptrs.h"
#include "algo/2d/SimpleStraightSkel.h"
#include "algo/3d/KernelWrapper.h"
#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "data/2d/Polygon.h"
#include "data/2d/skel/Arc.h"
#include "data/2d/skel/StraightSkeleton.h"
#include "data/2d/skel/SkelEdgeData.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SkelEdgeData.h"
#include <list>

namespace algo { namespace _3d {

WeightVertexSplitter::WeightVertexSplitter() {
    type_ = AbstractVertexSplitter::WEIGHT_VERTEX_SPLITTER;
    controller_ = ControllerSPtr();
}

WeightVertexSplitter::WeightVertexSplitter(ControllerSPtr controller) {
    type_ = AbstractVertexSplitter::WEIGHT_VERTEX_SPLITTER;
    controller_ = controller;
}

WeightVertexSplitter::~WeightVertexSplitter() {
    controller_.reset();
}

WeightVertexSplitterSPtr WeightVertexSplitter::create() {
    WeightVertexSplitterSPtr result =
            WeightVertexSplitterSPtr(new WeightVertexSplitter());
    return result;
}

WeightVertexSplitterSPtr WeightVertexSplitter::create(ControllerSPtr controller) {
    WeightVertexSplitterSPtr result =
            WeightVertexSplitterSPtr(new WeightVertexSplitter(controller));
    return result;
}

Plane3SPtr WeightVertexSplitter::createIntersectionPlane(VertexSPtr vertex) {
    Plane3SPtr result;
    Vector3SPtr normal = KernelFactory::createVector3(0.0, 0.0, 0.0);
    std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
    while (it_e != vertex->edges().end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge(edge_wptr);
            Vector3SPtr dir = KernelFactory::createVector3(edge->line());;
            if (edge->getVertexSrc() != vertex) {
                dir = KernelFactory::createVector3((*dir) * (-1));
            }
            normal = KernelFactory::createVector3((*normal) + (*dir));
        }
    }
    normal = KernelFactory::createVector3((*normal) / vertex->degree());
    normal = KernelWrapper::normalize(normal);
    Point3SPtr p = KernelFactory::createPoint3((*(vertex->getPoint())) + (*normal));
    normal = KernelFactory::createVector3((*normal) * (-1));

//    Vector3SPtr rotation = KernelFactory::createVector3(0.5, 0.0, 0.0);
//    normal = KernelFactory::createVector3((*normal) + (*rotation));
//    normal = KernelWrapper::normalize(normal);

    result = KernelFactory::createPlane3(p, normal);
    return result;
}

FacetSPtr WeightVertexSplitter::intersectPolyhedron(VertexSPtr vertex, Plane3SPtr plane) {
    FacetSPtr result = Facet::create();
    result->setPlane(plane);

    FacetSPtr polygon;
    EdgeSPtr edge_begin = vertex->firstEdge();
    EdgeSPtr edge_next = edge_begin->next(vertex);
    if (edge_begin->getFacetL() == edge_next->getFacetL() ||
            edge_begin->getFacetL() == edge_next->getFacetR()) {
        polygon = edge_begin->getFacetL();
    } else {
        polygon = edge_begin->getFacetR();
    }

    VertexSPtr vertex_begin;
    VertexSPtr vertex_src;
    VertexSPtr vertex_dst;
    EdgeSPtr edge_src;
    EdgeSPtr edge_dst;
    while (edge_src != edge_begin) {
        if (!edge_src) {
            edge_src = edge_begin;
            Point3SPtr p_src = KernelWrapper::intersection(
                    plane, edge_src->line());
            vertex_src = Vertex::create(p_src);
            result->addVertex(vertex_src);
            vertex_begin = vertex_src;
        }
        edge_dst = edge_src->next(vertex);
        if (edge_dst == edge_begin) {
            vertex_dst = vertex_begin;
        } else {
            Point3SPtr p_dst = KernelWrapper::intersection(
                    plane,edge_dst->line());
            vertex_dst = Vertex::create(p_dst);
            result->addVertex(vertex_dst);
        }

        EdgeSPtr edge = Edge::create(vertex_src, vertex_dst);
        edge->setFacetL(result);
        skel::SkelEdgeDataSPtr data = skel::SkelEdgeData::create(edge);
        data->setFacetOrigin(polygon);
        result->addEdge(edge);

        vertex_src = vertex_dst;
        edge_src = edge_dst;
        polygon = polygon->next(vertex);
    }

    return result;
}

data::_2d::PolygonSPtr WeightVertexSplitter::toWeighted2d(FacetSPtr polygon) {
    data::_2d::PolygonSPtr result = polygon->toPolygon();
    std::list<EdgeSPtr>::iterator it_e = polygon->edges().begin();
    std::list<data::_2d::EdgeSPtr>::iterator it_e2 = result->edges().begin();
    while (it_e != polygon->edges().end() &&
            it_e2 != result->edges().end()) {
        EdgeSPtr edge = *it_e++;
        data::_2d::EdgeSPtr edge2 = *it_e2++;
        skel::SkelEdgeDataSPtr data = std::dynamic_pointer_cast<skel::SkelEdgeData>(
                edge->getData());
        double angle = KernelWrapper::angle(
                polygon->plane(), data->getFacetOrigin()->plane());
        double speed = 1.0/sin(angle);
        data::_2d::skel::SkelEdgeDataSPtr data2 =
                data::_2d::skel::SkelEdgeData::create(edge2);
        data2->setSpeed(speed);
    }
    return result;
}

PolyhedronSPtr WeightVertexSplitter::splitVertex(VertexSPtr vertex) {
    PolyhedronSPtr polyhedron = vertex->getPolyhedron();
    if (vertex->degree() <= 3) {
        return polyhedron;
    }
    vertex->getData()->setHighlight(true);
    if (controller_) {
        controller_->wait();
    }
    Plane3SPtr plane = createIntersectionPlane(vertex);
    FacetSPtr polygon = intersectPolyhedron(vertex, plane);
    data::_2d::PolygonSPtr polygon_2d = toWeighted2d(polygon);
    DEBUG_VAR(polygon_2d->toString());
    algo::_2d::SimpleStraightSkelSPtr algoskel_2d =
            algo::_2d::SimpleStraightSkel::create(polygon_2d, controller_);
    data::_2d::skel::StraightSkeletonSPtr skel_2d = algoskel_2d->getResult();
    if (controller_) {
        controller_->setDispPolygon(polygon_2d);
        controller_->setDispSkel2d(skel_2d);
    }
    algoskel_2d->run();
    if (controller_) {
        controller_->haltSkip();
        controller_->wait();
    }

    WriteLock l(polyhedron->mutex());

    std::list<VertexSPtr> vertices_tosplit;
    vertices_tosplit.push_back(vertex);

    std::list<data::_2d::skel::ArcSPtr>::iterator it_a = skel_2d->arcs().begin();
    while (it_a != skel_2d->arcs().end()) {
        data::_2d::skel::ArcSPtr arc = *it_a++;
        data::_2d::EdgeSPtr edge_l = arc->getEdgeLeft();
        data::_2d::EdgeSPtr edge_r = arc->getEdgeRight();
        FacetSPtr polygon_l;
        FacetSPtr polygon_r;
        std::list<EdgeSPtr>::iterator it_e = polygon->edges().begin();
        std::list<data::_2d::EdgeSPtr>::iterator it_e2 = polygon_2d->edges().begin();
        while (it_e != polygon->edges().end() &&
                it_e2 != polygon_2d->edges().end()) {
            EdgeSPtr edge = *it_e++;
            data::_2d::EdgeSPtr edge_2d = *it_e2++;
            if (edge_2d == edge_l) {
                skel::SkelEdgeDataSPtr data_l = std::dynamic_pointer_cast<skel::SkelEdgeData>(
                        edge->getData());
                polygon_l = data_l->getFacetOrigin();
            }
            if (edge_2d == edge_r) {
                skel::SkelEdgeDataSPtr data_r = std::dynamic_pointer_cast<skel::SkelEdgeData>(
                        edge->getData());
                polygon_r = data_r->getFacetOrigin();
            }
        }

        if (polygon_l->findEdge(polygon_r)) {
            continue;
        }

        std::list<VertexSPtr>::iterator it_v = vertices_tosplit.begin();
        while (it_v != vertices_tosplit.end()) {
            std::list<VertexSPtr>::iterator it_current = it_v;
            VertexSPtr vertex_tosplit = *it_v++;
            if (vertex_tosplit->containsFacet(polygon_l) &&
                    vertex_tosplit->containsFacet(polygon_r)) {
                vertices_tosplit.erase(it_current);
                VertexSPtr vertex_dst = vertex_tosplit->split(polygon_l, polygon_r);
                if (vertex_tosplit->hasData()) {
                    skel::SkelVertexDataSPtr data_tosplit = std::dynamic_pointer_cast<skel::SkelVertexData>(
                            vertex_tosplit->getData());
                    skel::SkelVertexDataSPtr data_dst = skel::SkelVertexData::create(vertex_dst);
                    data_dst->setNode(data_tosplit->getNode());
                }
                if (vertex_tosplit->degree() > 3) {
                    vertices_tosplit.push_back(vertex_tosplit);
                }
                if (vertex_dst->degree() > 3) {
                    vertices_tosplit.push_back(vertex_dst);
                }
                break;
            }
        }
    }

    vertex->getData()->setHighlight(false);
    return polyhedron;
}

} }
