/**
 * @file   algo/3d/RotSimpleSphericalSkel.cpp
 * @author Gernot Walzl
 * @date   2012-12-28
 */

#include "algo/3d/RotSimpleSphericalSkel.h"

#include "debug.h"
#include "typedefs_thread.h"
#include "algo/Controller.h"
#include "algo/3d/KernelWrapper.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/SphericalPolygon.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/skel/SphericalSkeleton.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/CircularArc.h"
#include "data/3d/skel/SphericalSkelVertexData.h"
#include "data/3d/skel/SphericalSkelEdgeData.h"
#include "data/3d/skel/SphericalAbstractEvent.h"
#include "data/3d/skel/SphericalConstOffsetEvent.h"
#include "data/3d/skel/SphericalEdgeEvent.h"
#include "data/3d/skel/SphericalSplitEvent.h"
#include "data/3d/skel/SphericalTriangleEvent.h"
#include "util/Configuration.h"
#include <limits>

namespace algo { namespace _3d {

RotSimpleSphericalSkel::RotSimpleSphericalSkel(SphericalPolygonSPtr polygon) {
    type_ = AbstractSimpleSphericalSkel::ROT_SIMPLE_SPHERICAL_SKEL;
    polygon_ = polygon;
    controller_ = ControllerSPtr();
    skel_result_ = SphericalSkeleton::create(polygon->getSphere());
}

RotSimpleSphericalSkel::RotSimpleSphericalSkel(SphericalPolygonSPtr polygon, ControllerSPtr controller) {
    polygon_ = polygon;
    controller_ = controller;
    skel_result_ = SphericalSkeleton::create(polygon->getSphere());
}

RotSimpleSphericalSkel::~RotSimpleSphericalSkel() {
    // intentionally does nothing
}

RotSimpleSphericalSkelSPtr RotSimpleSphericalSkel::create(SphericalPolygonSPtr polygon) {
    return RotSimpleSphericalSkelSPtr(new RotSimpleSphericalSkel(polygon));
}

RotSimpleSphericalSkelSPtr RotSimpleSphericalSkel::create(SphericalPolygonSPtr polygon, ControllerSPtr controller) {
    return RotSimpleSphericalSkelSPtr(new RotSimpleSphericalSkel(polygon, controller));
}

void RotSimpleSphericalSkel::run() {
    if (controller_) {
        controller_->wait();
    }
    DEBUG_PRINT("== Rotational Spherical Skeleton started ==");
    SphericalPolygonSPtr polygon = polygon_;
    unsigned int i = 0;
    if (init(polygon)) {
        if (controller_) {
            controller_->wait();
        }
        double offset = 0.0;
        double offset_prev = 0.0;
        SphericalAbstractEventSPtr event = nextEvent(polygon, offset);
        while (event) {
            DEBUG_VAL("-- Next Event: " << event->toString() << " --");
            if (controller_) {
                controller_->wait();
            }
            offset = event->getOffset();
            polygon = shiftEdges(polygon, offset - offset_prev);
            if (event->getType() == SphericalAbstractEvent::CONST_OFFSET_EVENT) {
                event->setPolygonResult(polygon);
                skel_result_->addEvent(event);
            } else if (event->getType() == SphericalAbstractEvent::EDGE_EVENT) {
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
            event = nextEvent(polygon, offset);
            offset_prev = offset;
        }
        DEBUG_PRINT("== Spherical Skeleton finished ==");
        DEBUG_VAR(skel_result_->toString());
    }
}


bool RotSimpleSphericalSkel::isReflex(CircularVertexSPtr vertex) {
    bool result = false;
    CircularEdgeSPtr edge_in = vertex->getEdgeIn();
    CircularEdgeSPtr edge_out = vertex->getEdgeOut();
    Sphere3SPtr sphere = vertex->getPolygon()->getSphere();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    Vector3SPtr dir = KernelFactory::createVector3(*vertex->getPoint() - *p_center);
    Vector3SPtr normal_in = KernelFactory::createVector3(edge_in->supportingPlane());
    Vector3SPtr v_trans = KernelWrapper::cross(normal_in, dir);
    Point3SPtr p_trans = KernelFactory::createPoint3(*p_center + *v_trans);
    if (KernelWrapper::side(edge_out->supportingPlane(), p_trans) < 0) {
        result = true;
    }
    return result;
}


bool RotSimpleSphericalSkel::init(SphericalPolygonSPtr polygon) {
    WriteLock l(polygon->mutex());
    bool result = true;
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
            CircularNodeSPtr node = CircularNode::create(vertex->getPoint());
            data->setNode(node);
            CircularArcSPtr arc = createArc(vertex);
            skel_result_->addNode(node);
            skel_result_->addArc(arc);
        } else {
            result = false;
        }
    }
    initSpeeds(polygon);
    return result;
}


void RotSimpleSphericalSkel::initSpeeds(SphericalPolygonSPtr polygon) {
    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        CircularEdgeSPtr edge_in = vertex->getEdgeIn();
        CircularEdgeSPtr edge_out = vertex->getEdgeOut();
        Plane3SPtr plane_in = edge_in->supportingPlane();
        Plane3SPtr plane_out = edge_out->supportingPlane();
        double angle = M_PI - KernelWrapper::angle(plane_in, plane_out);
        if (isReflex(vertex)) {
            angle = 2.0*M_PI - angle;
        }
        double speed = 1.0/sin(angle/2.0);
        SphericalSkelVertexDataSPtr data;
        if (vertex->hasData()) {
            data = std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
        } else {
            data = SphericalSkelVertexData::create(vertex);
        }
        data->setSpeed(speed);
    }
}


void RotSimpleSphericalSkel::updateSpeeds(SphericalPolygonSPtr polygon_prev) {
    Sphere3SPtr sphere = polygon_prev->getSphere();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    std::list<CircularVertexSPtr>::iterator it_v = polygon_prev->vertices().begin();
    while (it_v != polygon_prev->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        Vector3SPtr dir = KernelFactory::createVector3(
                *(vertex->getPoint()) - *p_center);
        SphericalSkelVertexDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
        CircularVertexSPtr vertex_rotated = data->getOffsetVertex();
        Vector3SPtr dir_rotated = KernelFactory::createVector3(
                *(vertex_rotated->getPoint()) - *p_center);
        double angle = KernelWrapper::angle(dir, dir_rotated);
        data->setSpeed(angle);
    }
}


double RotSimpleSphericalSkel::approxOffsetTo(CircularVertexSPtr vertex, Point3SPtr point) {
    double result = 0.0;
    Sphere3SPtr sphere = vertex->getPolygon()->getSphere();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    Vector3SPtr dir_point = KernelFactory::createVector3(
            *point - *p_center);
    Vector3SPtr dir_vertex = KernelFactory::createVector3(
            *(vertex->getPoint()) - *p_center);
    double angle = KernelWrapper::angle(dir_point, dir_vertex);
    SphericalSkelVertexDataSPtr data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
    result = -angle / data->getSpeed();
    return result;
}


SphericalEdgeEventSPtr RotSimpleSphericalSkel::nextEdgeEvent(SphericalPolygonSPtr polygon, double offset) {
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
        double offset_src = approxOffsetTo(vertex_src, point);
        double offset_dst = approxOffsetTo(vertex_dst, point);
        double offset_current = (offset_src + offset_dst)/2.0;
        if (offset_max < offset_current && offset_current <= 0.0) {
            CircularNodeSPtr node;
            if (!result) {
                node = CircularNode::create(point);
                result = SphericalEdgeEvent::create();
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(offset_current + offset);
            node->setPoint(point);
            result->setEdge(edge);
            SphericalSkelVertexDataSPtr data_src =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(
                    edge->getVertexSrc()->getData());
            SphericalSkelVertexDataSPtr data_dst =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(
                    edge->getVertexDst()->getData());
            node->addArc(data_src->getArc());
            node->addArc(data_dst->getArc());
            offset_max = offset_current;
        }
    }
    return result;
}

SphericalSplitEventSPtr RotSimpleSphericalSkel::nextSplitEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalSplitEventSPtr result;
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        if (!isReflex(vertex)) {
            // no split event
            continue;
        }
        std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
        while (it_e != polygon->edges().end()) {
            CircularEdgeSPtr edge = *it_e++;
            if (edge == vertex->getEdgeIn() ||
                    edge == vertex->getEdgeOut()) {
                continue;
            }
            if (edge->next()->getVertexDst() == vertex ||
                    edge->prev()->getVertexSrc() == vertex) {
                // sticking event
                continue;
            }
            Point3SPtr point = crashAt(vertex, edge);
            if (!point) {
                continue;
            }
            double offset_current = approxOffsetTo(vertex, point);
            if (offset_max < offset_current && offset_current <= 0.0) {
                CircularNodeSPtr node;
                if (!result) {
                    node = CircularNode::create(point);
                    result = SphericalSplitEvent::create();
                    result->setNode(node);
                }
                node = result->getNode();
                node->clear();
                node->setOffset(offset_current + offset);
                node->setPoint(point);
                result->setVertex(vertex);
                result->setEdge(edge);
                SphericalSkelVertexDataSPtr data =
                        std::dynamic_pointer_cast<SphericalSkelVertexData>(
                        vertex->getData());
                node->addArc(data->getArc());
                offset_max = offset_current;
            }
        }
    }
    return result;
}

SphericalTriangleEventSPtr RotSimpleSphericalSkel::nextTriangleEvent(SphericalPolygonSPtr polygon, double offset) {
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
        CircularVertexSPtr vertex_dst = edge->getVertexDst();
        double offset_src = approxOffsetTo(vertex_src, point);
        double offset_dst = approxOffsetTo(vertex_dst, point);
        double offset_current = (offset_src + offset_dst)/2.0;
        if (offset_max < offset_current && offset_current <= 0.0) {
            CircularNodeSPtr node;
            if (!result) {
                node = CircularNode::create(point);
                result = SphericalTriangleEvent::create();
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(offset_current + offset);
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


SphericalAbstractEventSPtr RotSimpleSphericalSkel::nextEvent(SphericalPolygonSPtr polygon, double offset) {
    SphericalAbstractEventSPtr result = SphericalAbstractEventSPtr();
    if (!polygon) {
        return result;
    }
    if (polygon->edges().size() == 0) {
        return result;
    }
    SphericalAbstractEventSPtr events[4];
    for (unsigned int i = 0; i < 4; i++) {
        events[i] = SphericalAbstractEventSPtr();
    }
    double const_offset = util::Configuration::getInstance()->getDouble(
            "algo_3d_RotSimpleSphericalSkel", "const_offset");
    if (const_offset != 0.0) {
        events[0] = SphericalConstOffsetEvent::create(offset + const_offset);
    }
    events[1] = nextEdgeEvent(polygon, offset);
    events[2] = nextSplitEvent(polygon, offset);
    events[3] = nextTriangleEvent(polygon, offset);
    for (unsigned int i = 0; i < 4; i++) {
        if (events[i]) {
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


CircularEdgeSPtr RotSimpleSphericalSkel::findLongestEdge(std::list<CircularEdgeSPtr> edges) {
    CircularEdgeSPtr result;
    double angle_max = 0.0;
    Sphere3SPtr sphere = polygon_->getSphere();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    std::list<CircularEdgeSPtr>::iterator it_e = edges.begin();
    while (it_e != edges.end()) {
        CircularEdgeSPtr edge = *it_e++;
        Vector3SPtr dir_src = KernelFactory::createVector3(
                *(edge->getVertexSrc()->getPoint()) - *p_center);
        Vector3SPtr dir_dst = KernelFactory::createVector3(
                *(edge->getVertexDst()->getPoint()) - *p_center);
        double angle = KernelWrapper::angle(dir_src, dir_dst);
        if (angle > angle_max) {
            result = edge;
            angle_max = angle;
        }
    }
    return result;
}


double RotSimpleSphericalSkel::distanceOffset(CircularEdgeSPtr edge_begin, double offset, double speed_dst) {
    double result = 0.0;
    SphericalPolygonSPtr polygon = edge_begin->getPolygon();
    Sphere3SPtr sphere = polygon->getSphere();
    double radius = polygon->getRadius();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    Point3SPtr p_begin_rot;
    Point3SPtr p_src_rot;
    Point3SPtr p_dst_rot;
    Plane3SPtr plane_edge_rot;
    Plane3SPtr plane_edge_rot_prev;
    CircularEdgeSPtr edge;
    while (edge != edge_begin) {
        if (!edge) {
            edge = edge_begin;

            CircularVertexSPtr vertex_src = edge->getVertexSrc();
            SphericalSkelVertexDataSPtr data_src =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_src->getData());
            CircularArcSPtr arc_src = data_src->getArc();
            Vector3SPtr axis_src = KernelWrapper::normalize(
                    KernelFactory::createVector3(arc_src->getSupportingPlane()));
            Vector3SPtr dir_src = KernelWrapper::normalize(
                    KernelFactory::createVector3(*(vertex_src->getPoint()) - *p_center));
            Vector3SPtr dir_src_rot = KernelWrapper::rotateVector(
                    dir_src, axis_src, -offset * data_src->getSpeed());
            p_src_rot = KernelFactory::createPoint3(
                    *p_center + ((*dir_src_rot)*radius));
            p_begin_rot = p_src_rot;

            CircularVertexSPtr vertex_dst = edge->getVertexDst();
            SphericalSkelVertexDataSPtr data_dst =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_dst->getData());
            CircularArcSPtr arc_dst = data_dst->getArc();
            Vector3SPtr axis_dst = KernelWrapper::normalize(
                    KernelFactory::createVector3(arc_dst->getSupportingPlane()));
            Vector3SPtr dir_dst = KernelWrapper::normalize(
                    KernelFactory::createVector3(*(vertex_dst->getPoint()) - *p_center));
            Vector3SPtr dir_dst_rot = KernelWrapper::rotateVector(
                    dir_dst, axis_dst, -offset * speed_dst);
            p_dst_rot = KernelFactory::createPoint3(
                        *p_center + ((*dir_dst_rot)*radius));

            plane_edge_rot = KernelFactory::createPlane3(p_center, p_src_rot, p_dst_rot);
        } else {
            CircularVertexSPtr vertex_src = edge->getVertexSrc();
            SphericalSkelVertexDataSPtr data_src =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_src->getData());
            Plane3SPtr plane_src_arc = data_src->getArc()->getSupportingPlane();
            double angle = M_PI - KernelWrapper::angle(plane_src_arc, plane_edge_rot_prev);
            Line3SPtr axis = KernelWrapper::intersection(plane_src_arc, plane_edge_rot_prev);
            plane_edge_rot = KernelWrapper::rotatePlane(plane_src_arc, axis, angle);

            CircularVertexSPtr vertex_dst = edge->getVertexDst();
            SphericalSkelVertexDataSPtr data_dst =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_dst->getData());
            CircularArcSPtr arc_dst = data_dst->getArc();
            Plane3SPtr plane_dst_arc = arc_dst->getSupportingPlane();
            Line3SPtr line_dst_rot = KernelWrapper::intersection(plane_dst_arc, plane_edge_rot);
            Vector3SPtr dir_dst_rot = KernelWrapper::normalize(KernelFactory::createVector3(line_dst_rot));
            p_dst_rot = KernelFactory::createPoint3(
                    *p_center + ((*dir_dst_rot)*radius));
        }
        p_src_rot = p_dst_rot;
        plane_edge_rot_prev = plane_edge_rot;
        edge = edge->next();
    }
    result = KernelWrapper::distance(p_dst_rot, p_begin_rot);
    return result;
}


double RotSimpleSphericalSkel::findMinDistance(CircularEdgeSPtr edge_begin, double offset) {
    double result = 0.0;
    CircularVertexSPtr vertex_dst = edge_begin->getVertexDst();
    SphericalSkelVertexDataSPtr data_dst =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_dst->getData());
    double speed = data_dst->getSpeed();
    double dist = 0.0;
    double dist_min = std::numeric_limits<double>::max();
    for (double speed_factor = 0.25; speed_factor <= 4.0; speed_factor += 0.001) {
        dist = distanceOffset(edge_begin, offset, speed * speed_factor);
        if (dist < dist_min) {
            result = speed * speed_factor;
            dist_min = dist;
        }
    }
    return result;
}


SphericalPolygonSPtr RotSimpleSphericalSkel::shiftEdges(SphericalPolygonSPtr polygon, double offset) {
    Sphere3SPtr sphere = polygon->getSphere();
    SphericalPolygonSPtr result = SphericalPolygon::create(sphere);
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    double radius = polygon->getRadius();

    std::list<CircularEdgeSPtr> edges_torotate;
    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        edges_torotate.push_back(edge);
    }

    while (!edges_torotate.empty()) {
        CircularVertexSPtr vertex_begin_rot;
        CircularVertexSPtr vertex_src_rot;
        CircularVertexSPtr vertex_dst_rot;
        CircularEdgeSPtr edge_begin = findLongestEdge(edges_torotate);
        double speed_dst = findMinDistance(edge_begin, offset);
        CircularVertexSPtr vertex_dst = edge_begin->getVertexDst();
        SphericalSkelVertexDataSPtr data_dst =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_dst->getData());
        DEBUG_VAR(data_dst->getSpeed());
        data_dst->setSpeed(speed_dst);
        DEBUG_VAR(data_dst->getSpeed());
        CircularEdgeSPtr edge;
        while (edge != edge_begin) {
            if (!edge) {
                edge = edge_begin;

                CircularVertexSPtr vertex_src = edge->getVertexSrc();
                SphericalSkelVertexDataSPtr data_src =
                        std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_src->getData());
                data_src->setHighlight(true);
                CircularArcSPtr arc_src = data_src->getArc();
                Vector3SPtr axis_src = KernelWrapper::normalize(
                        KernelFactory::createVector3(arc_src->getSupportingPlane()));
                Vector3SPtr dir_src = KernelWrapper::normalize(
                        KernelFactory::createVector3(*(vertex_src->getPoint()) - *p_center));
                Vector3SPtr dir_src_rot = KernelWrapper::rotateVector(
                        dir_src, axis_src, -offset * data_src->getSpeed());
                Point3SPtr p_src_rot = KernelFactory::createPoint3(
                        *p_center + ((*dir_src_rot)*radius));
                vertex_src_rot = CircularVertex::create(p_src_rot);
                SphericalSkelVertexDataSPtr data_src_rot =
                        SphericalSkelVertexData::create(vertex_src_rot);
                data_src_rot->setArc(arc_src);
                data_src_rot->setSpeed(data_src->getSpeed());
                data_src->setOffsetVertex(vertex_src_rot);
                result->addVertex(vertex_src_rot);
                vertex_begin_rot = vertex_src_rot;

                CircularVertexSPtr vertex_dst = edge->getVertexDst();
                SphericalSkelVertexDataSPtr data_dst =
                        std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_dst->getData());
                CircularArcSPtr arc_dst = data_dst->getArc();
                Vector3SPtr axis_dst = KernelWrapper::normalize(
                        KernelFactory::createVector3(arc_dst->getSupportingPlane()));
                Vector3SPtr dir_dst = KernelWrapper::normalize(
                        KernelFactory::createVector3(*(vertex_dst->getPoint()) - *p_center));
                Vector3SPtr dir_dst_rot = KernelWrapper::rotateVector(
                        dir_dst, axis_dst, -offset * data_dst->getSpeed());
                Point3SPtr p_dst_rot = KernelFactory::createPoint3(
                        *p_center + ((*dir_dst_rot)*radius));
                vertex_dst_rot = CircularVertex::create(p_dst_rot);
                SphericalSkelVertexDataSPtr data_dst_rot =
                        SphericalSkelVertexData::create(vertex_dst_rot);
                data_dst_rot->setArc(arc_dst);
                data_dst_rot->setSpeed(data_dst->getSpeed());
                data_dst->setOffsetVertex(vertex_dst_rot);
                result->addVertex(vertex_dst_rot);
            } else {
                CircularVertexSPtr vertex_src = edge->getVertexSrc();
                SphericalSkelVertexDataSPtr data_src =
                        std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_src->getData());
                assert(vertex_src_rot == data_src->getOffsetVertex());
                Plane3SPtr plane_src_arc = data_src->getArc()->getSupportingPlane();
                Plane3SPtr plane_src_in = vertex_src_rot->getEdgeIn()->supportingPlane();
                double angle = M_PI - KernelWrapper::angle(plane_src_arc, plane_src_in);
                Line3SPtr axis = KernelWrapper::intersection(plane_src_arc, plane_src_in);
                Plane3SPtr plane_edge_rot = KernelWrapper::rotatePlane(plane_src_arc, axis, angle);

                CircularVertexSPtr vertex_dst = edge->getVertexDst();
                SphericalSkelVertexDataSPtr data_dst =
                        std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_dst->getData());
                CircularArcSPtr arc_dst = data_dst->getArc();
                Plane3SPtr plane_dst_arc = arc_dst->getSupportingPlane();
                Line3SPtr line_dst_rot = KernelWrapper::intersection(plane_dst_arc, plane_edge_rot);
                Vector3SPtr dir_dst_rot = KernelWrapper::normalize(KernelFactory::createVector3(line_dst_rot));
                Point3SPtr p_dst_rot = KernelFactory::createPoint3(
                        *p_center + ((*dir_dst_rot)*radius));

                if (edge->next() == edge_begin) {
                    vertex_dst_rot = vertex_begin_rot;
                    DEBUG_VAR(KernelWrapper::distance(p_dst_rot, vertex_begin_rot->getPoint()));
                } else {
                    vertex_dst_rot = CircularVertex::create(p_dst_rot);
                    SphericalSkelVertexDataSPtr data_dst_rot =
                            SphericalSkelVertexData::create(vertex_dst_rot);
                    data_dst_rot->setArc(arc_dst);
                    data_dst_rot->setSpeed(data_dst->getSpeed());
                    data_dst->setOffsetVertex(vertex_dst_rot);
                    result->addVertex(vertex_dst_rot);
                }
            }

            CircularEdgeSPtr edge_rot = CircularEdge::create(vertex_src_rot, vertex_dst_rot);
            SphericalSkelEdgeDataSPtr data;
            if (edge->hasData()) {
                data = std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge->getData());
            } else {
                data = SphericalSkelEdgeData::create(edge);
            }
            data->setOffsetEdge(edge_rot);
            result->addEdge(edge_rot);

            edges_torotate.remove(edge);
            vertex_src_rot = vertex_dst_rot;
            edge = edge->next();
        }
    }

    updateSpeeds(polygon);
    return result;
}


void RotSimpleSphericalSkel::checkAngles(SphericalPolygonSPtr polygon) {
    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        SphericalSkelVertexDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
        Plane3SPtr plane_arc = data->getArc()->getSupportingPlane();
        Plane3SPtr plane_in = vertex->getEdgeIn()->supportingPlane();
        Plane3SPtr plane_out = vertex->getEdgeOut()->supportingPlane();
        double angle_in =  M_PI - KernelWrapper::angle(plane_arc, plane_in);
        double angle_out = KernelWrapper::angle(plane_arc, plane_out);
        DEBUG_VAR(angle_in - angle_out);
    }
}


SphericalPolygonSPtr RotSimpleSphericalSkel::shiftEdges2(SphericalPolygonSPtr polygon, double offset) {
    Sphere3SPtr sphere = polygon->getSphere();
    SphericalPolygonSPtr result = SphericalPolygon::create(sphere);
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    double radius = polygon->getRadius();

    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        SphericalSkelVertexDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
        CircularArcSPtr arc = data->getArc();
        Vector3SPtr axis = KernelWrapper::normalize(
                KernelFactory::createVector3(arc->getSupportingPlane()));
        Vector3SPtr dir = KernelWrapper::normalize(
                KernelFactory::createVector3(*(vertex->getPoint()) - *p_center));
        Vector3SPtr dir_rotated = KernelWrapper::rotateVector(
                dir, axis, -offset * data->getSpeed());
        Point3SPtr p_rotated = KernelFactory::createPoint3(
                *p_center + ((*dir_rotated)*radius));
        CircularVertexSPtr vertex_rotated = CircularVertex::create(p_rotated);
        SphericalSkelVertexDataSPtr data_rotated =
                SphericalSkelVertexData::create(vertex_rotated);
        data_rotated->setArc(arc);
        data_rotated->setSpeed(data->getSpeed());
        data->setOffsetVertex(vertex_rotated);
        result->addVertex(vertex_rotated);
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
        CircularVertexSPtr vertex_src_rotated = data_src->getOffsetVertex();
        CircularVertexSPtr vertex_dst_rotated = data_dst->getOffsetVertex();
        CircularEdgeSPtr edge_rotated = CircularEdge::create(
                vertex_src_rotated, vertex_dst_rotated);
        data->setOffsetEdge(edge_rotated);
        result->addEdge(edge_rotated);
    }

    initSpeeds(result);
    return result;
}


void RotSimpleSphericalSkel::handleEdgeEvent(SphericalEdgeEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularNodeSPtr node = event->getNode();
    appendEventNode(node);

    // remove edge and link adjacent edges
    SphericalSkelEdgeDataSPtr edge_data =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(
            event->getEdge()->getData());
    CircularEdgeSPtr edge_toremove = edge_data->getOffsetEdge();
    CircularVertexSPtr vertex = edge_toremove->getVertexSrc();
    CircularVertexSPtr vertex_dst = edge_toremove->getVertexDst();
    SphericalSkelVertexDataSPtr vertex_data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(
            vertex->getData());
    polygon->removeEdge(edge_toremove);
    vertex->setEdgeOut(vertex_dst->getEdgeOut());
    vertex_dst->getEdgeOut()->setVertexSrc(vertex);
    vertex_dst->setEdgeOut(CircularEdgeSPtr());
    polygon->removeVertex(vertex_dst);
    vertex->setPoint(node->getPoint());

    // create arc
    vertex_data->setNode(node);
    CircularArcSPtr arc = createArc(vertex);

    // save and return result
    skel_result_->addArc(arc);
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void RotSimpleSphericalSkel::handleSplitEvent(SphericalSplitEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularNodeSPtr node = event->getNode();
    appendEventNode(node);

    // split edge
    SphericalSkelVertexDataSPtr vertex_data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(
            event->getVertex()->getData());
    CircularVertexSPtr vertex = vertex_data->getOffsetVertex();
    vertex->setPoint(node->getPoint());
    SphericalSkelEdgeDataSPtr edge_data =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(
            event->getEdge()->getData());
    CircularEdgeSPtr edge = edge_data->getOffsetEdge();
    CircularVertexSPtr vertex_dst = edge->getVertexDst();
    CircularEdgeSPtr edge_l = vertex->getEdgeIn();

    edge->setVertexDst(vertex);
    vertex->setEdgeIn(edge);
    SphericalSkelVertexDataSPtr vertex_data_r =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(
            vertex->getData());
    CircularVertexSPtr vertex_l = CircularVertex::create(node->getPoint());
    SphericalSkelVertexDataSPtr vertex_data_l =
            SphericalSkelVertexData::create(vertex_l);
    polygon->addVertex(vertex_l);
    edge_l->setVertexDst(vertex_l);
    vertex_l->setEdgeIn(edge_l);
    CircularEdgeSPtr edge_2 = CircularEdge::create(vertex_l, vertex_dst);
    SphericalSkelEdgeDataSPtr edge_data_2 = SphericalSkelEdgeData::create(edge_2);
    edge_data_2->setSpeed(edge_data->getSpeed());
    edge_data_2->setRotationAxis(edge_data->getRotationAxis());
    edge_data_2->setFacetOrigin(edge_data->getFacetOrigin());
    polygon->addEdge(edge_2);

    // create arcs
    vertex_data_l->setNode(node);
    vertex_data_r->setNode(node);
    CircularArcSPtr arc_l = createArc(vertex_l);
    CircularArcSPtr arc_r = createArc(vertex);

    // save and return result
    skel_result_->addArc(arc_l);
    skel_result_->addArc(arc_r);
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void RotSimpleSphericalSkel::handleTriangleEvent(SphericalTriangleEventSPtr event, SphericalPolygonSPtr polygon) {
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
        CircularEdgeSPtr offset_edge = data->getOffsetEdge();
        polygon->removeEdge(offset_edge);
    }

    // remove vertices
    CircularVertexSPtr vertices[3];
    event->getVertices(vertices);
    for (unsigned int i = 0; i < 3; i++) {
        SphericalSkelVertexDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                vertices[i]->getData());
        CircularVertexSPtr offset_vertex = data->getOffsetVertex();
        polygon->removeVertex(offset_vertex);
    }

    // save and return result
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}


} }
