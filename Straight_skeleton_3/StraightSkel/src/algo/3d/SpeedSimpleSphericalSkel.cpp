/**
 * @file   algo/3d/SpeedSimpleSphericalSkel.cpp
 * @author Gernot Walzl
 * @date   2013-09-04
 */

#include "algo/3d/SpeedSimpleSphericalSkel.h"

#include "debug.h"
#include "algo/Controller.h"
#include "algo/3d/KernelWrapper.h"
#include "data/3d/SphericalPolygon.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/Facet.h"
#include "data/3d/skel/SphericalSkeleton.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/CircularArc.h"
#include "data/3d/skel/SphericalSkelVertexData.h"
#include "data/3d/skel/SphericalSkelEdgeData.h"
#include "data/3d/skel/SphericalAbstractEvent.h"
#include "data/3d/skel/SphericalEdgeEvent.h"
#include "data/3d/skel/SphericalSplitEvent.h"
#include "data/3d/skel/SphericalTriangleEvent.h"
#include "util/Configuration.h"
#include <list>

namespace algo { namespace _3d {

SpeedSimpleSphericalSkel::SpeedSimpleSphericalSkel(SphericalPolygonSPtr polygon) {
    type_ = AbstractSimpleSphericalSkel::SPEED_SIMPLE_SPHERICAL_SKEL;
    polygon_ = polygon;
    controller_ = ControllerSPtr();
    skel_result_ = SphericalSkeleton::create(polygon->getSphere());
}

SpeedSimpleSphericalSkel::SpeedSimpleSphericalSkel(SphericalPolygonSPtr polygon, ControllerSPtr controller) {
    type_ = AbstractSimpleSphericalSkel::SPEED_SIMPLE_SPHERICAL_SKEL;
    polygon_ = polygon;
    controller_ = controller;
    skel_result_ = SphericalSkeleton::create(polygon->getSphere());
}

SpeedSimpleSphericalSkel::~SpeedSimpleSphericalSkel() {
    // intentionally does nothing
}

SpeedSimpleSphericalSkelSPtr SpeedSimpleSphericalSkel::create(SphericalPolygonSPtr polygon) {
    return SpeedSimpleSphericalSkelSPtr(new SpeedSimpleSphericalSkel(polygon));
}

SpeedSimpleSphericalSkelSPtr SpeedSimpleSphericalSkel::create(SphericalPolygonSPtr polygon, ControllerSPtr controller) {
    return SpeedSimpleSphericalSkelSPtr(new SpeedSimpleSphericalSkel(polygon, controller));
}

void SpeedSimpleSphericalSkel::run() {
    if (controller_) {
        controller_->wait();
    }
    DEBUG_PRINT("== Speed Spherical Skeleton started ==");
    SphericalPolygonSPtr polygon = polygon_;
    unsigned int i = 0;
    if (init(polygon)) {
        if (controller_) {
            controller_->wait();
        }
        SphericalAbstractEventSPtr event = nextEvent(polygon);
        while (event) {
            DEBUG_VAL("-- Next Event: " << event->toString() << " --");
            if (controller_) {
                controller_->wait();
            }
            polygon = copyPolygon(polygon);
            if (event->getType() == SphericalAbstractEvent::EDGE_EVENT) {
                handleEdgeEvent(std::dynamic_pointer_cast<SphericalEdgeEvent>(event), polygon);
            } else if (event->getType() == SphericalAbstractEvent::SPLIT_EVENT) {
                handleSplitEvent(std::dynamic_pointer_cast<SphericalSplitEvent>(event), polygon);
            } else if (event->getType() == SphericalAbstractEvent::TRIANGLE_EVENT) {
                handleTriangleEvent(std::dynamic_pointer_cast<SphericalTriangleEvent>(event), polygon);
            }
            assert(polygon->isConsistent());
            assert(skel_result_->isConsistent());
            DEBUG_PRINT("-- Finished handling Event --");
            i++;
            DEBUG_VAR(i);
            if (controller_) {
                controller_->wait();
            }
            event = nextEvent(polygon);
        }
        DEBUG_PRINT("== Spherical Skeleton finished ==");
        DEBUG_VAR(skel_result_->toString());
    }
}


bool SpeedSimpleSphericalSkel::isReflex(CircularVertexSPtr vertex) {
    bool result = false;
    SphericalSkelVertexDataSPtr data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
    result = data->isReflex();
    return result;
}

bool SpeedSimpleSphericalSkel::init(SphericalPolygonSPtr polygon) {
    WriteLock l(polygon->mutex());
    bool result = true;
    Sphere3SPtr sphere = polygon->getSphere();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        CircularEdgeSPtr edge_in = vertex->getEdgeIn();
        CircularEdgeSPtr edge_out = vertex->getEdgeOut();
        if (edge_in && edge_out) {
            SphericalSkelVertexDataSPtr data;
            if (vertex->hasData()) {
                data = std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
            } else {
                data = SphericalSkelVertexData::create(vertex);
            }
            Vector3SPtr dir = KernelFactory::createVector3(*vertex->getPoint() - *p_center);
            Vector3SPtr normal_in = KernelFactory::createVector3(edge_in->supportingPlane());
            Vector3SPtr v_trans = KernelWrapper::cross(normal_in, dir);
            Point3SPtr p_trans = KernelFactory::createPoint3(*p_center + *v_trans);
            if (KernelWrapper::side(edge_out->supportingPlane(), p_trans) < 0) {
                data->setReflex(true);
            } else {
                data->setReflex(false);
            }
            CircularNodeSPtr node = CircularNode::create(vertex->getPoint());
            data->setNode(node);
            CircularArcSPtr arc = createArc(vertex);
            data->setSpeed(speed(vertex));
            skel_result_->addNode(node);
            skel_result_->addArc(arc);
        } else {
            result = false;
        }
    }
    //initSpeeds(polygon);
    return result;
}


void SpeedSimpleSphericalSkel::initSpeeds(SphericalPolygonSPtr polygon) {
    Sphere3SPtr sphere = polygon->getSphere();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    CircularEdgeSPtr first;
    CircularEdgeSPtr edge = polygon->edges().front();
    while (edge != first) {
        CircularVertexSPtr vertex_src = edge->getVertexSrc();
        CircularVertexSPtr vertex_dst = edge->getVertexDst();
        SphericalSkelVertexDataSPtr data_src =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                vertex_src->getData());
        SphericalSkelVertexDataSPtr data_dst =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                vertex_dst->getData());
        if (!first) {
            data_src->setSpeed(1.0);
            first = edge;
        }
        Point3SPtr p_vanish = vanishesAt(edge);
        Vector3SPtr dir_vanish = KernelFactory::createVector3(*p_vanish - *p_center);
        Vector3SPtr dir_src = KernelFactory::createVector3(*(vertex_src->getPoint()) - *p_center);
        Vector3SPtr dir_dst = KernelFactory::createVector3(*(vertex_dst->getPoint()) - *p_center);
        double angle_src = KernelWrapper::angle(dir_src, dir_vanish);
        double angle_dst = KernelWrapper::angle(dir_dst, dir_vanish);
        data_dst->setSpeed(data_src->getSpeed()*angle_dst/angle_src);
        edge = edge->next();
    }
}


double SpeedSimpleSphericalSkel::speed(CircularVertexSPtr vertex) {
    double result = 0.0;
    CircularEdgeSPtr edge_in = getEdgeOrigin(vertex->getEdgeIn());
    CircularEdgeSPtr edge_out = getEdgeOrigin(vertex->getEdgeOut());
    double angle = M_PI - KernelWrapper::angle(
            edge_in->supportingPlane(), edge_out->supportingPlane());
    result = 1.0/sin(angle/2.0);
    return result;
}


Point3SPtr SpeedSimpleSphericalSkel::vanishesAt(CircularEdgeSPtr edge) {
    Point3SPtr result;
    SphericalPolygonSPtr polygon = edge->getPolygon();
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    double radius = polygon->getRadius();
    CircularEdgeSPtr edge_origin = getEdgeOrigin(edge);
    SphericalSkelEdgeDataSPtr data_origin =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(
            edge_origin->getData());
    Plane3SPtr plane_origin = data_origin->getFacetOrigin()->plane();
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
        Point3SPtr p_vanish = KernelFactory::createPoint3((*p_center)
                + ((*dir) * radius));
        if (KernelWrapper::side(plane_origin, p_vanish) > 0) {
            p_vanish = KernelFactory::createPoint3((*p_center)
                    - ((*dir) * radius));
        }
        result = p_vanish;
    }
    return result;
}

Point3SPtr SpeedSimpleSphericalSkel::crashAt(CircularVertexSPtr vertex, CircularEdgeSPtr edge) {
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
        Line3SPtr line = KernelWrapper::intersection(bisector_out, bisector_in);
        if (line) {
            Vector3SPtr dir = KernelWrapper::normalize(
                    KernelFactory::createVector3(line));
            Point3SPtr p_crash = KernelFactory::createPoint3((*p_center)
                    + ((*dir) * radius));
            if (KernelWrapper::side(plane_origin, p_crash) > 0) {
                p_crash = KernelFactory::createPoint3((*p_center)
                        - ((*dir) * radius));
            }
        }
    }
    if (p_crash) {
        // check if p_crash is inside bounds of edge
        CircularVertexSPtr vertex_src = edge->getVertexSrc();
        CircularVertexSPtr vertex_dst = edge->getVertexDst();
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
            }
        }
    }
    return result;
}


double SpeedSimpleSphericalSkel::offsetTo(CircularVertexSPtr vertex, Point3SPtr point) {
    double result = 0.0;
    Point3SPtr p_center = KernelFactory::createPoint3(vertex->getPolygon()->getSphere());
    Vector3SPtr dir_vertex = KernelFactory::createVector3(*(vertex->getPoint()) - *p_center);
    Vector3SPtr dir_point = KernelFactory::createVector3(*point - *p_center);
    double angle = KernelWrapper::angle(dir_vertex, dir_point);
    SphericalSkelVertexDataSPtr data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
    double speed = data->getSpeed();
    result = - angle / speed;
    return result;
}


SphericalEdgeEventSPtr SpeedSimpleSphericalSkel::nextEdgeEvent(SphericalPolygonSPtr polygon) {
    ReadLock l(polygon->mutex());
    SphericalEdgeEventSPtr result;
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        if (isTriangle(edge)) {
            // triangle event
            continue;
        }
        Point3SPtr point = vanishesAt(edge);
        if (!point) {
            continue;
        }
        CircularVertexSPtr vertex_src = edge->getVertexSrc();
        CircularVertexSPtr vertex_dst = edge->getVertexDst();
        double offset_src = offsetTo(vertex_src, point);
        double offset_dst = offsetTo(vertex_dst, point);
        double accuracy_edge_event = offset_dst - offset_src;
        DEBUG_VAR(accuracy_edge_event);
        SphericalSkelVertexDataSPtr data_src =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                vertex_src->getData());
        SphericalSkelVertexDataSPtr data_dst =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                vertex_dst->getData());
        CircularNodeSPtr node_src = data_src->getNode();
        double offset_current = node_src->getOffset() + offset_src;
        if (offset_max < offset_current) {
            CircularNodeSPtr node;
            if (!result) {
                node = CircularNode::create(point);
                result = SphericalEdgeEvent::create();
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(offset_current);
            node->setPoint(point);
            result->setEdge(edge);
            node->addArc(data_src->getArc());
            node->addArc(data_dst->getArc());
            offset_max = offset_current;
        }
    }
    return result;
}

SphericalSplitEventSPtr SpeedSimpleSphericalSkel::nextSplitEvent(SphericalPolygonSPtr polygon) {
    ReadLock l(polygon->mutex());
    SphericalSplitEventSPtr result;
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
        while (it_e != polygon->edges().end()) {
            CircularEdgeSPtr edge = *it_e++;
            if (edge == vertex->getEdgeIn() ||
                    edge == vertex->getEdgeOut()) {
                continue;
            }
            if (vertex->prev() == edge->getVertexDst() ||
                    vertex->next() == edge->getVertexSrc()) {
                continue;
            }
            Point3SPtr point = crashAt(vertex, edge);
            if (!point) {
                continue;
            }
            double offset_vertex = offsetTo(vertex, point);
            SphericalSkelVertexDataSPtr vertex_data =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(
                    vertex->getData());
            CircularNodeSPtr node_src = vertex_data->getNode();
            double offset_current = node_src->getOffset() + offset_vertex;
            if (offset_max < offset_current) {
                CircularNodeSPtr node;
                if (!result) {
                    node = CircularNode::create(point);
                    result = SphericalSplitEvent::create();
                    result->setNode(node);
                }
                node = result->getNode();
                node->clear();
                node->setOffset(offset_current);
                node->setPoint(point);
                result->setVertex(vertex);
                result->setEdge(edge);
                node->addArc(vertex_data->getArc());
                offset_max = offset_current;
            }
        }
    }
    return result;
}

SphericalTriangleEventSPtr SpeedSimpleSphericalSkel::nextTriangleEvent(SphericalPolygonSPtr polygon) {
    ReadLock l(polygon->mutex());
    SphericalTriangleEventSPtr result;
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        if (!isTriangle(edge)) {
            // edge event
            continue;
        }
        Point3SPtr point = vanishesAt(edge);
        if (!point) {
            continue;
        }
        CircularVertexSPtr vertex_src = edge->getVertexSrc();
        double offset_src = offsetTo(vertex_src, point);
        SphericalSkelVertexDataSPtr data_src =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                vertex_src->getData());
        CircularNodeSPtr node = data_src->getNode();
        double offset_current = node->getOffset() + offset_src;
        if (offset_max < offset_current) {
            CircularNodeSPtr node;
            if (!result) {
                node = CircularNode::create(point);
                result = SphericalTriangleEvent::create();
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(offset_current);
            node->setPoint(point);
            result->setEdgeBegin(edge);
            for (unsigned int i = 0; i < 3; i++) {
                SphericalSkelVertexDataSPtr data =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(
                    edge->getVertexSrc()->getData());
                node->addArc(data->getArc());
                edge = edge->next();
            }
            offset_max = offset_current;
        }
    }
    return result;
}


SphericalAbstractEventSPtr SpeedSimpleSphericalSkel::nextEvent(SphericalPolygonSPtr polygon) {
    SphericalAbstractEventSPtr result = SphericalAbstractEventSPtr();
    if (!polygon) {
        return result;
    }
    if (polygon->edges().size() == 0) {
        return result;
    }
    SphericalAbstractEventSPtr events[3];
    for (unsigned int i = 0; i < 3; i++) {
        events[i] = SphericalAbstractEventSPtr();
    }
    events[0] = nextEdgeEvent(polygon);
    events[1] = nextSplitEvent(polygon);
    events[2] = nextTriangleEvent(polygon);
    for (unsigned int i = 0; i < 3; i++) {
        if (events[i]) {
            DEBUG_VAR(events[i]->toString());
            if (!result) {
                result = events[i];
            } else if (result->getOffset() <= events[i]->getOffset()) {
                result = events[i];
            }
        }
    }
    result->setHighlight(true);
    return result;
}


SphericalPolygonSPtr SpeedSimpleSphericalSkel::copyPolygon(SphericalPolygonSPtr polygon) {
    Sphere3SPtr sphere = polygon->getSphere();
    SphericalPolygonSPtr result = SphericalPolygon::create(sphere);

    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        CircularVertexSPtr vertex_c = vertex->clone();
        SphericalSkelVertexDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
        SphericalSkelVertexDataSPtr data_c =
                SphericalSkelVertexData::create(vertex_c);
        data_c->setNode(data->getNode());
        data_c->setArc(data->getArc());
        data_c->setReflex(data->isReflex());
        data_c->setSpeed(data->getSpeed());
        data->setOffsetVertex(vertex_c);
        result->addVertex(vertex_c);
    }

    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        SphericalSkelEdgeDataSPtr data;
        if (edge->hasData()) {
            data = std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge->getData());
        } else {
            data = SphericalSkelEdgeData::create(edge);
        }
        SphericalSkelVertexDataSPtr data_src =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                edge->getVertexSrc()->getData());
        SphericalSkelVertexDataSPtr data_dst =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                edge->getVertexDst()->getData());
        CircularVertexSPtr vertex_src_c = data_src->getOffsetVertex();
        CircularVertexSPtr vertex_dst_c = data_dst->getOffsetVertex();
        if (vertex_src_c && vertex_dst_c) {
            CircularEdgeSPtr edge_c = CircularEdge::create(
                    vertex_src_c, vertex_dst_c);
            edge_c->setSupportingPlane(edge->getSupportingPlane());
            data->setOffsetEdge(edge_c);
            result->addEdge(edge_c);
        } else {
            DEBUG_VAR(edge->toString());
        }
    }

    return result;
}


void SpeedSimpleSphericalSkel::handleEdgeEvent(SphericalEdgeEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularNodeSPtr node = event->getNode();
    appendEventNode(node);

    // remove edge and link adjacent edges
    SphericalSkelEdgeDataSPtr edge_data =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(
            event->getEdge()->getData());
    CircularEdgeSPtr edge = edge_data->getOffsetEdge();

    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    CircularVertexSPtr vertex = edge->getVertexSrc();
    CircularVertexSPtr vertex_dst = edge->getVertexDst();

    SphericalSkelVertexDataSPtr vertex_data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(
            vertex->getData());
    polygon->removeEdge(edge);
    vertex->setEdgeOut(vertex_dst->getEdgeOut());
    vertex_dst->getEdgeOut()->setVertexSrc(vertex);
    vertex_dst->setEdgeOut(CircularEdgeSPtr());
    polygon->removeVertex(vertex_dst);
    vertex->setPoint(node->getPoint());

    // create arc
    vertex_data->setReflex(false);
    vertex_data->setNode(node);
    CircularArcSPtr arc = createArc(vertex);
    vertex_data->setSpeed(speed(vertex));

    CircularEdgeSPtr edge_in = vertex->getEdgeIn();
    CircularEdgeSPtr edge_out = vertex->getEdgeOut();
    Plane3SPtr plane_in = KernelFactory::createPlane3(p_center,
            edge_in->getVertexSrc()->getPoint(), edge_in->getVertexDst()->getPoint());
    edge_in->setSupportingPlane(plane_in);
    Plane3SPtr plane_out = KernelFactory::createPlane3(p_center,
            edge_out->getVertexSrc()->getPoint(), edge_out->getVertexDst()->getPoint());
    edge_out->setSupportingPlane(plane_out);

    // save and return result
    skel_result_->addArc(arc);
    skel_result_->addEvent(event);
    event->setPolygonResult(polygon);
}

void SpeedSimpleSphericalSkel::handleSplitEvent(SphericalSplitEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularNodeSPtr node = event->getNode();
    appendEventNode(node);

    // split edge
    // TODO

    // save and return result
    skel_result_->addEvent(event);
    event->setPolygonResult(polygon);
}

void SpeedSimpleSphericalSkel::handleTriangleEvent(SphericalTriangleEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularNodeSPtr node = event->getNode();
    appendEventNode(node);

    // remove edges
    CircularEdgeSPtr edges[3];
    event->getEdges(edges);
    for (unsigned int i = 0; i < 3; i++) {
        SphericalSkelEdgeDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelEdgeData>(
                edges[i]->getData());
        CircularEdgeSPtr edge = data->getOffsetEdge();
        polygon->removeEdge(edge);
    }

    // remove vertices
    CircularVertexSPtr vertices[3];
    event->getVertices(vertices);
    for (unsigned int i = 0; i < 3; i++) {
        SphericalSkelVertexDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                vertices[i]->getData());
        CircularVertexSPtr vertex = data->getOffsetVertex();
        polygon->removeVertex(vertex);
    }

    // save and return result
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

} }
