/**
 * @file   algo/3d/ProjSimpleSphericalSkel.cpp
 * @author Gernot Walzl
 * @date   2012-12-03
 */

#include "algo/3d/ProjSimpleSphericalSkel.h"

#include "debug.h"
#include "typedefs_thread.h"
#include "algo/Controller.h"
#include "algo/3d/KernelWrapper.h"
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
#include <list>
#include <limits>

namespace algo { namespace _3d {

ProjSimpleSphericalSkel::ProjSimpleSphericalSkel(SphericalPolygonSPtr polygon) {
    type_ = AbstractSimpleSphericalSkel::PROJ_SIMPLE_SPHERICAL_SKEL;
    polygon_ = polygon;
    controller_ = ControllerSPtr();
    skel_result_ = SphericalSkeleton::create(polygon->getSphere());
}

ProjSimpleSphericalSkel::ProjSimpleSphericalSkel(SphericalPolygonSPtr polygon, ControllerSPtr controller) {
    polygon_ = polygon;
    controller_ = controller;
    skel_result_ = SphericalSkeleton::create(polygon->getSphere());
}

ProjSimpleSphericalSkel::~ProjSimpleSphericalSkel() {
    // intentionally does nothing
}

ProjSimpleSphericalSkelSPtr ProjSimpleSphericalSkel::create(SphericalPolygonSPtr polygon) {
    return ProjSimpleSphericalSkelSPtr(new ProjSimpleSphericalSkel(polygon));
}

ProjSimpleSphericalSkelSPtr ProjSimpleSphericalSkel::create(SphericalPolygonSPtr polygon, ControllerSPtr controller) {
    return ProjSimpleSphericalSkelSPtr(new ProjSimpleSphericalSkel(polygon, controller));
}

void ProjSimpleSphericalSkel::run() {
    if (controller_) {
        controller_->wait();
    }
    DEBUG_PRINT("== Projective Spherical Skeleton started ==");
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


bool ProjSimpleSphericalSkel::isReflex(CircularVertexSPtr vertex) {
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


void ProjSimpleSphericalSkel::initRotationAxes(SphericalPolygonSPtr polygon) {
    Sphere3SPtr sphere = polygon->getSphere();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    Vector3SPtr normal = KernelFactory::createVector3(0.0, 0.0, 1.0);  // TODO
    Plane3SPtr plane_rot_axes = KernelFactory::createPlane3(p_center, normal);
    skel_result_->setRotAxesPlane(plane_rot_axes);
    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        Point3SPtr p_src = edge->getVertexSrc()->getPoint();
        Point3SPtr p_dst = edge->getVertexDst()->getPoint();
        int side_src = KernelWrapper::side(plane_rot_axes, p_src);
        int side_dst = KernelWrapper::side(plane_rot_axes, p_dst);
        Plane3SPtr plane_edge = edge->supportingPlane();
        Line3SPtr rotation_axis;
        if (side_src < 0 && side_dst < 0) {
            rotation_axis = KernelWrapper::intersection(
                    plane_edge, plane_rot_axes);
        } else if (side_src > 0 && side_dst > 0) {
            rotation_axis = KernelWrapper::intersection(
                    plane_rot_axes, plane_edge);
        }
        SphericalSkelEdgeDataSPtr data;
        if (edge->hasData()) {
            data = std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge->getData());
        } else {
            data = SphericalSkelEdgeData::create(edge);
        }
        data->setRotationAxis(rotation_axis);
    }
}

void ProjSimpleSphericalSkel::initConstSpeeds(SphericalPolygonSPtr polygon) {
    Plane3SPtr plane_rot_axes = skel_result_->getRotAxesPlane();
    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        CircularEdgeSPtr edge_origin = getEdgeOrigin(edge);
        Plane3SPtr plane_edge = edge_origin->supportingPlane();
        double angle = KernelWrapper::angle(plane_rot_axes, plane_edge);
        double speed = 1.0 / sin(angle);
        SphericalSkelEdgeDataSPtr data;
        if (edge->hasData()) {
            data = std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge->getData());
        } else {
            data = SphericalSkelEdgeData::create(edge);
        }
        data->setSpeed(speed);
    }
}

bool ProjSimpleSphericalSkel::init(SphericalPolygonSPtr polygon) {
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
    initRotationAxes(polygon);
    initConstSpeeds(polygon);
    return result;
}


bool ProjSimpleSphericalSkel::hasValidRotationAxis(CircularEdgeSPtr edge) {
    bool result = false;
    if (!edge) {
        return result;
    }
    Plane3SPtr plane_rot_axes = skel_result_->getRotAxesPlane();
    int side_src = KernelWrapper::side(plane_rot_axes,
            edge->getVertexSrc()->getPoint());
    int side_dst = KernelWrapper::side(plane_rot_axes,
            edge->getVertexDst()->getPoint());
    result = (side_src == side_dst);
    return result;
}


double ProjSimpleSphericalSkel::angleTo(CircularEdgeSPtr edge) {
    double result = 0.0;
    Plane3SPtr plane_rot_axes = skel_result_->getRotAxesPlane();
    Plane3SPtr plane_edge = edge->supportingPlane();
    result = KernelWrapper::angle(plane_rot_axes, plane_edge);
    return result;
}


double ProjSimpleSphericalSkel::offsetTo(CircularEdgeSPtr edge, Point3SPtr point) {
    double result = -std::numeric_limits<double>::infinity();
    if (!hasValidRotationAxis(edge)) {
        return result;
    }
    double angle = angleTo(edge);
    Sphere3SPtr sphere = edge->getPolygon()->getSphere();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    SphericalSkelEdgeDataSPtr data =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(
            edge->getData());
    Line3SPtr rotation_axis = data->getRotationAxis();
    Plane3SPtr plane_edge = edge->supportingPlane();
    Vector3SPtr dir = KernelFactory::createVector3(rotation_axis);
    Point3SPtr p_center_moved = KernelFactory::createPoint3(*p_center + *dir);
    Plane3SPtr plane_point = KernelFactory::createPlane3(p_center, p_center_moved, point);
    double delta_angle = KernelWrapper::angle(plane_edge, plane_point);
    if (KernelWrapper::side(plane_edge, point) < 0) {
        delta_angle *= -1.0;
    }
    if ((angle - delta_angle) > 0.0) {
        double offset = sin(delta_angle) / (sin(angle) * sin(angle - delta_angle));
        if (offset <= 0.0) {
            result = offset / data->getSpeed();
        }
    }
    return result;
}


SphericalEdgeEventSPtr ProjSimpleSphericalSkel::nextEdgeEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalEdgeEventSPtr result;
    double offset_max = -std::numeric_limits<double>::infinity();
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
        double offset_current = offsetTo(edge, point);
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

SphericalSplitEventSPtr ProjSimpleSphericalSkel::nextSplitEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalSplitEventSPtr result;
    double offset_max = -std::numeric_limits<double>::infinity();
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
            double offset_current = offsetTo(edge, point);
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

SphericalTriangleEventSPtr ProjSimpleSphericalSkel::nextTriangleEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalTriangleEventSPtr result;
    double offset_max = -std::numeric_limits<double>::infinity();
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
        double offset_current = offsetTo(edge, point);
        if (offset_max <= offset_current && offset_current <= 0.0) {
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

SphericalAbstractEventSPtr ProjSimpleSphericalSkel::nextEvent(SphericalPolygonSPtr polygon, double offset) {
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
            "algo_3d_ProjSimpleSphericalSkel", "const_offset");
    if (const_offset != 0.0) {
        events[0] = SphericalConstOffsetEvent::create(offset + const_offset);
    }
    events[1] = nextEdgeEvent(polygon, offset);
    events[2] = nextSplitEvent(polygon, offset);
    events[3] = nextTriangleEvent(polygon, offset);
    if (!events[1] && !events[2] && events[3]) {
        events[0] = SphericalAbstractEventSPtr();
    }
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


SphericalPolygonSPtr ProjSimpleSphericalSkel::shiftEdges(SphericalPolygonSPtr polygon, double offset) {
    Sphere3SPtr sphere = polygon->getSphere();
    SphericalPolygonSPtr result = SphericalPolygon::create(sphere);
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    double radius = polygon->getRadius();

    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        CircularEdgeSPtr edge_in = vertex->getEdgeIn();
        CircularEdgeSPtr edge_out = vertex->getEdgeOut();
        bool has_valid_rot_axis_in = hasValidRotationAxis(edge_in);
        bool has_valid_rot_axis_out = hasValidRotationAxis(edge_out);
        SphericalSkelVertexDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
        CircularArcSPtr arc = data->getArc();
        CircularVertexSPtr vertex_rotated;
        Point3SPtr p_in;
        Point3SPtr p_out;
        if (has_valid_rot_axis_in) {
            Plane3SPtr plane_arc = arc->getSupportingPlane();
            double angle = angleTo(edge_in);
            Plane3SPtr plane_edge = edge_in->supportingPlane();
            SphericalSkelEdgeDataSPtr data_edge =
                    std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge_in->getData());
            Line3SPtr rotation_axis = data_edge->getRotationAxis();
            if (angle > M_PI/2.0) {
                angle = M_PI - angle;
            }
            double diff = (1.0/tan(angle)) + (offset * data_edge->getSpeed());
            double delta_angle = 0.0;
            if (diff == 0.0) {
                delta_angle = -(M_PI/2.0)+angle;
            } else if (diff < 0.0) {
                delta_angle = atan(diff) - (M_PI/2.0)+angle;
            } else {
                delta_angle = -atan(1.0/diff) + angle;
            }
            Plane3SPtr plane_rotated = KernelWrapper::rotatePlane(
                    plane_edge, rotation_axis, delta_angle);
            Line3SPtr line_rotated = KernelWrapper::intersection(
                    plane_arc, plane_rotated);
            Vector3SPtr dir_rotated = KernelWrapper::normalize(
                    KernelFactory::createVector3(line_rotated));
            p_in = KernelFactory::createPoint3(
                    *p_center + ((*dir_rotated) * radius));
        }
        if (has_valid_rot_axis_out) {
            Plane3SPtr plane_arc = arc->getSupportingPlane();
            double angle = angleTo(edge_out);
            Plane3SPtr plane_edge = edge_out->supportingPlane();
            SphericalSkelEdgeDataSPtr data_edge =
                    std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge_out->getData());
            Line3SPtr rotation_axis = data_edge->getRotationAxis();
            if (angle > M_PI/2.0) {
                angle = M_PI - angle;
            }
            double diff = (1.0/tan(angle)) + (offset * data_edge->getSpeed());
            double delta_angle = 0.0;
            if (diff == 0.0) {
                delta_angle = -(M_PI/2.0)+angle;
            } else if (diff < 0.0) {
                delta_angle = atan(diff) - (M_PI/2.0)+angle;
            } else {
                delta_angle = -atan(1.0/diff) + angle;
            }
            Plane3SPtr plane_rotated = KernelWrapper::rotatePlane(
                    plane_edge, rotation_axis, delta_angle);
            Line3SPtr line_rotated = KernelWrapper::intersection(
                    plane_arc, plane_rotated);
            Vector3SPtr dir_rotated = KernelWrapper::normalize(
                    KernelFactory::createVector3(line_rotated));
            p_out = KernelFactory::createPoint3(
                    *p_center + ((*dir_rotated) * radius));
        }
        if (has_valid_rot_axis_in && has_valid_rot_axis_out) {
            DEBUG_VAR(KernelWrapper::distance(p_in, p_out));
            vertex_rotated = CircularVertex::create(p_in);
        }
        if (has_valid_rot_axis_in && !has_valid_rot_axis_out) {
            vertex_rotated = CircularVertex::create(p_in);
        }
        if (!has_valid_rot_axis_in && has_valid_rot_axis_out) {
            vertex_rotated = CircularVertex::create(p_out);
        }
        if (!has_valid_rot_axis_in && !has_valid_rot_axis_out) {
            DEBUG_VAL("Error: CircularVertex does not have an adjacent CircularEdge with a valid rotation axis.");
            DEBUG_VAR(vertex->toString());
            vertex_rotated = CircularVertex::create(vertex->getPoint());
            vertex_rotated->setPointValid(false);
        }
        SphericalSkelVertexDataSPtr data_rotated =
                SphericalSkelVertexData::create(vertex_rotated);
        data_rotated->setArc(arc);
        data->setOffsetVertex(vertex_rotated);
        result->addVertex(vertex_rotated);
    }

    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        SphericalSkelEdgeDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelEdgeData>(edge->getData());
        SphericalSkelVertexDataSPtr data_src =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                edge->getVertexSrc()->getData());
        SphericalSkelVertexDataSPtr data_dst =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(
                edge->getVertexDst()->getData());
        CircularVertexSPtr vertex_src_rotated = data_src->getOffsetVertex();
        CircularVertexSPtr vertex_dst_rotated = data_dst->getOffsetVertex();
        if (vertex_src_rotated && vertex_dst_rotated) {
            CircularEdgeSPtr edge_rotated = CircularEdge::create(
                    vertex_src_rotated, vertex_dst_rotated);
            SphericalSkelEdgeDataSPtr data_rotated =
                    SphericalSkelEdgeData::create(edge_rotated);
            data_rotated->setSpeed(data->getSpeed());
            data_rotated->setRotationAxis(data->getRotationAxis());
            data->setOffsetEdge(edge_rotated);
            result->addEdge(edge_rotated);
        } else {
            DEBUG_VAL("Error: CircularEdge has an adjacent CircularVertex that was not rotated.");
            DEBUG_VAR(edge->toString());
        }
    }

    return result;
}


void ProjSimpleSphericalSkel::handleEdgeEvent(SphericalEdgeEventSPtr event, SphericalPolygonSPtr polygon) {
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

void ProjSimpleSphericalSkel::handleSplitEvent(SphericalSplitEventSPtr event, SphericalPolygonSPtr polygon) {
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

void ProjSimpleSphericalSkel::handleTriangleEvent(SphericalTriangleEventSPtr event, SphericalPolygonSPtr polygon) {
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
