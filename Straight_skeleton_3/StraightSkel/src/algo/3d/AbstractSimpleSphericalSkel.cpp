/**
 * @file   algo/3d/AbstractSimpleSphericalSkel.cpp
 * @author Gernot Walzl
 * @date   2013-01-11
 */

#include "algo/3d/AbstractSimpleSphericalSkel.h"

#include "debug.h"
#include "algo/Controller.h"
#include "algo/3d/KernelWrapper.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/Facet.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/SphericalPolygon.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/CircularArc.h"
#include "data/3d/skel/SphericalSkeleton.h"
#include "data/3d/skel/SphericalSkelVertexData.h"
#include "data/3d/skel/SphericalSkelEdgeData.h"
#include <list>

namespace algo { namespace _3d {

AbstractSimpleSphericalSkel::AbstractSimpleSphericalSkel() {
    type_ = 0;
}

AbstractSimpleSphericalSkel::~AbstractSimpleSphericalSkel() {
    polygon_.reset();
    controller_.reset();
    skel_result_.reset();
}


int AbstractSimpleSphericalSkel::getType() const {
    return type_;
}


ThreadSPtr AbstractSimpleSphericalSkel::startThread() {
    return ThreadSPtr(new std::thread(
            std::bind(&AbstractSimpleSphericalSkel::run, this)));
}


CircularEdgeSPtr AbstractSimpleSphericalSkel::getEdgeOrigin(CircularEdgeSPtr edge) {
    CircularEdgeSPtr result = edge;
    CircularVertexSPtr vertex_dst = edge->getVertexDst();
    if (vertex_dst->hasData()) {
        SphericalSkelVertexDataSPtr data_dst =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(
            vertex_dst->getData());
        CircularArcSPtr arc_dst = data_dst->getArc();
        if (arc_dst) {
            result = arc_dst->getEdgeLeft();
        }
    }
    CircularVertexSPtr vertex_src = edge->getVertexSrc();
    if (vertex_src->hasData()) {
        SphericalSkelVertexDataSPtr data_src =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(
            vertex_src->getData());
        CircularArcSPtr arc_src = data_src->getArc();
        if (arc_src) {
            result = arc_src->getEdgeRight();
        }
    }
    return result;
}


CircularArcSPtr AbstractSimpleSphericalSkel::createArc(CircularVertexSPtr vertex) {
    CircularArcSPtr result;
    SphericalSkelVertexDataSPtr data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
    data->setArc(CircularArcSPtr());
    CircularNodeSPtr node = data->getNode();
    CircularEdgeSPtr edge_l = getEdgeOrigin(vertex->getEdgeIn());
    CircularEdgeSPtr edge_r = getEdgeOrigin(vertex->getEdgeOut());
    SphericalSkelEdgeDataSPtr data_l =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge_l->getData());
    SphericalSkelEdgeDataSPtr data_r =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge_r->getData());
    Plane3SPtr plane_l_origin = data_l->getFacetOrigin()->plane();
    Plane3SPtr plane_r_origin = data_r->getFacetOrigin()->plane();
    Plane3SPtr plane_bisector;
    if (isReflex(vertex)) {
        plane_bisector = KernelWrapper::bisector(
                plane_l_origin, KernelWrapper::opposite(plane_r_origin));
    } else {
        plane_bisector = KernelWrapper::bisector(
                KernelWrapper::opposite(plane_l_origin), plane_r_origin);
    }
    Sphere3SPtr sphere = vertex->getPolygon()->getSphere();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    Vector3SPtr dir_node = KernelFactory::createVector3(
            *(node->getPoint()) - *p_center);
    Vector3SPtr normal = KernelFactory::createVector3(plane_bisector);
    Vector3SPtr dir_arc = KernelWrapper::normalize(
            KernelWrapper::cross(normal, dir_node));
    result = CircularArc::create(node, dir_arc);
    result->setEdgeLeft(edge_l);
    result->setEdgeRight(edge_r);
    result->setSupportingPlane(plane_bisector);
    data->setArc(result);
    return result;
}


Point3SPtr AbstractSimpleSphericalSkel::vanishesAt(CircularEdgeSPtr edge) {
    Point3SPtr result;
    SphericalPolygonSPtr polygon = edge->getPolygon();
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    double radius = polygon->getRadius();
    SphericalSkelVertexDataSPtr data_src =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(
            edge->getVertexSrc()->getData());
    SphericalSkelVertexDataSPtr data_dst =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(
            edge->getVertexDst()->getData());
    CircularArcSPtr arc_src = data_src->getArc();
    CircularArcSPtr arc_dst = data_dst->getArc();
    Plane3SPtr plane_src = arc_src->getSupportingPlane();
    Plane3SPtr plane_dst = arc_dst->getSupportingPlane();
    Line3SPtr line = KernelWrapper::intersection(plane_dst, plane_src);
    if (line) {
        Vector3SPtr dir = KernelWrapper::normalize(
                KernelFactory::createVector3(line));
        result = KernelFactory::createPoint3((*p_center) + ((*dir) * radius));
    }
    return result;
}

Point3SPtr AbstractSimpleSphericalSkel::crashAt(CircularVertexSPtr vertex, CircularEdgeSPtr edge) {
    Point3SPtr result;
    CircularEdgeSPtr edge_in = vertex->getEdgeIn();
    CircularEdgeSPtr edge_out = vertex->getEdgeOut();
    if (edge_in == edge || edge_out == edge) {
        return result;
    }
    SphericalPolygonSPtr polygon = vertex->getPolygon();
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    double radius = polygon->getRadius();
    CircularEdgeSPtr edge_origin = getEdgeOrigin(edge);
    CircularEdgeSPtr edge_in_origin = getEdgeOrigin(edge_in);
    CircularEdgeSPtr edge_out_origin = getEdgeOrigin(edge_out);
    SphericalSkelEdgeDataSPtr data_origin =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(
            edge_origin->getData());
    SphericalSkelEdgeDataSPtr data_in_origin =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(
            edge_in_origin->getData());
    SphericalSkelEdgeDataSPtr data_out_origin =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(
            edge_out_origin->getData());
    Plane3SPtr plane_origin = data_origin->getFacetOrigin()->plane();
    Plane3SPtr plane_in_origin = data_in_origin->getFacetOrigin()->plane();
    Plane3SPtr plane_out_origin = data_out_origin->getFacetOrigin()->plane();
    Plane3SPtr bisector_in = KernelWrapper::bisector(
            KernelWrapper::opposite(plane_in_origin), plane_origin);
    Plane3SPtr bisector_out = KernelWrapper::bisector(
            KernelWrapper::opposite(plane_origin), plane_out_origin);
    Point3SPtr p_crash;
    if (bisector_in && bisector_out) {
        bool is_reflex = true;
//        if (vertex->hasData()) {
//            SphericalSkelVertexDataSPtr data =
//                    std::dynamic_pointer_cast<SphericalSkelVertexData>(
//                    vertex->getData());
//            is_reflex = data->isReflex();
//        }
        Line3SPtr line;
        if (is_reflex) {
            line = KernelWrapper::intersection(bisector_out, bisector_in);
        } else {
            line = KernelWrapper::intersection(bisector_in, bisector_out);
        }
        if (line) {
            Vector3SPtr dir = KernelWrapper::normalize(
                    KernelFactory::createVector3(line));
            p_crash = KernelFactory::createPoint3((*p_center) + ((*dir) * radius));
        }
    }
    if (p_crash) {
        // check if p_crash is inside bounds of edge
        CircularVertexSPtr vertex_src = edge->getVertexSrc();
        CircularVertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->isPointValid() && vertex_dst->isPointValid()) {
            SphericalSkelVertexDataSPtr data_src =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(
                    vertex_src->getData());
            SphericalSkelVertexDataSPtr data_dst =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(
                    vertex_dst->getData());
            if (data_src && data_dst) {
                CircularArcSPtr arc_src = data_src->getArc();
                CircularArcSPtr arc_dst = data_dst->getArc();
                Plane3SPtr plane_src = arc_src->getSupportingPlane();
                Plane3SPtr plane_dst = arc_dst->getSupportingPlane();
                if (KernelWrapper::side(plane_src, p_crash) > 0 &&
                        KernelWrapper::side(plane_dst, p_crash) < 0) {
                    result = p_crash;
                } else if (KernelWrapper::side(plane_src, p_crash) > 0 &&
                        vertex->prev() == vertex_dst) {
                    result = p_crash;
                } else if (vertex->next() == vertex_src &&
                        KernelWrapper::side(plane_dst, p_crash) < 0) {
                    result = p_crash;
                } else if (vertex->next() == vertex_src &&
                        vertex->prev() == vertex_dst) {
                    result = p_crash;
                }
            }
        } else {
            result = p_crash;
        }
    }
    return result;
}


bool AbstractSimpleSphericalSkel::isTriangle(CircularEdgeSPtr edge_begin) {
    bool result = false;
    CircularEdgeSPtr edge = edge_begin;
    for (unsigned int i = 0; i < 3; i++) {
        edge = edge->next();
        if (!edge) {
            break;
        }
    }
    result = (edge == edge_begin);
    return result;
}


void AbstractSimpleSphericalSkel::appendEventNode(CircularNodeSPtr node) {
    std::list<CircularArcWPtr>::iterator it_a = node->arcs().begin();
    while (it_a != node->arcs().end()) {
        CircularArcWPtr arc_wptr = *it_a++;
        if (!arc_wptr.expired()) {
            CircularArcSPtr arc = CircularArcSPtr(arc_wptr);
            arc->setNodeDst(node);
            arc->setNodeDstListIt(
                    std::find(node->arcs().begin(), node->arcs().end(), arc_wptr));
        }
    }
    skel_result_->addNode(node);
}


SphericalSkeletonSPtr AbstractSimpleSphericalSkel::getResult() const {
    return this->skel_result_;
}

} }

