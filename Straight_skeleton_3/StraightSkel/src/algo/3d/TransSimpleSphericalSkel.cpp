/**
 * @file   algo/3d/TransSimpleSphericalSkel.cpp
 * @author Gernot Walzl
 * @date   2013-01-11
 */

#include "algo/3d/TransSimpleSphericalSkel.h"

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
#include "data/3d/skel/SphericalDblEdgeEvent.h"
#include "data/3d/skel/SphericalLeaveEvent.h"
#include "data/3d/skel/SphericalReturnEvent.h"
#include "data/3d/skel/SphericalDblLeaveEvent.h"
#include "data/3d/skel/SphericalDblReturnEvent.h"
#include "data/3d/skel/SphericalVertexEvent.h"
#include "data/3d/skel/SphericalEdgeMergeEvent.h"
#include "data/3d/skel/SphericalInversionEvent.h"
#include "util/Configuration.h"
#include <list>

namespace algo { namespace _3d {

TransSimpleSphericalSkel::TransSimpleSphericalSkel(SphericalPolygonSPtr polygon) {
    type_ = AbstractSimpleSphericalSkel::TRANS_SIMPLE_SPHERICAL_SKEL;
    polygon_ = polygon;
    controller_ = ControllerSPtr();
    skel_result_ = SphericalSkeleton::create(polygon->getSphere());
    epsilon_ = 0.005;
}

TransSimpleSphericalSkel::TransSimpleSphericalSkel(SphericalPolygonSPtr polygon, ControllerSPtr controller) {
    type_ = AbstractSimpleSphericalSkel::TRANS_SIMPLE_SPHERICAL_SKEL;
    polygon_ = polygon;
    controller_ = controller;
    skel_result_ = SphericalSkeleton::create(polygon->getSphere());
    epsilon_ = 0.005;
}

TransSimpleSphericalSkel::~TransSimpleSphericalSkel() {
    // intentionally does nothing
}

TransSimpleSphericalSkelSPtr TransSimpleSphericalSkel::create(SphericalPolygonSPtr polygon) {
    return TransSimpleSphericalSkelSPtr(new TransSimpleSphericalSkel(polygon));
}

TransSimpleSphericalSkelSPtr TransSimpleSphericalSkel::create(SphericalPolygonSPtr polygon, ControllerSPtr controller) {
    return TransSimpleSphericalSkelSPtr(new TransSimpleSphericalSkel(polygon, controller));
}

void TransSimpleSphericalSkel::run() {
    if (controller_) {
        controller_->wait();
    }
    DEBUG_PRINT("== Translational Spherical Skeleton started ==");
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
            } else if (event->getType() == SphericalAbstractEvent::DBL_EDGE_EVENT) {
                handleDblEdgeEvent(std::dynamic_pointer_cast<SphericalDblEdgeEvent>(event), polygon);
            } else if (event->getType() == SphericalAbstractEvent::LEAVE_EVENT) {
                handleLeaveEvent(std::dynamic_pointer_cast<SphericalLeaveEvent>(event), polygon);
            } else if (event->getType() == SphericalAbstractEvent::RETURN_EVENT) {
                handleReturnEvent(std::dynamic_pointer_cast<SphericalReturnEvent>(event), polygon);
            } else if (event->getType() == SphericalAbstractEvent::DBL_LEAVE_EVENT) {
                handleDblLeaveEvent(std::dynamic_pointer_cast<SphericalDblLeaveEvent>(event), polygon);
            } else if (event->getType() == SphericalAbstractEvent::DBL_RETURN_EVENT) {
                handleDblReturnEvent(std::dynamic_pointer_cast<SphericalDblReturnEvent>(event), polygon);
            } else if (event->getType() == SphericalAbstractEvent::VERTEX_EVENT) {
                handleVertexEvent(std::dynamic_pointer_cast<SphericalVertexEvent>(event), polygon);
            } else if (event->getType() == SphericalAbstractEvent::EDGE_MERGE_EVENT) {
                handleEdgeMergeEvent(std::dynamic_pointer_cast<SphericalEdgeMergeEvent>(event), polygon);
            } else if (event->getType() == SphericalAbstractEvent::INVERSION_EVENT) {
                handleInversionEvent(std::dynamic_pointer_cast<SphericalInversionEvent>(event), polygon);
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


bool TransSimpleSphericalSkel::isReflex(CircularVertexSPtr vertex) {
    bool result = false;
    SphericalSkelVertexDataSPtr data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
    result = data->isReflex();
    return result;
}


CircularArcSPtr TransSimpleSphericalSkel::createArc(CircularVertexSPtr vertex) {
    CircularArcSPtr result;
    SphericalSkelVertexDataSPtr data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
    data->setArc(CircularArcSPtr());
    CircularNodeSPtr node = data->getNode();
    CircularEdgeSPtr edge_l = vertex->getEdgeIn();
    CircularEdgeSPtr edge_r = vertex->getEdgeOut();
    Plane3SPtr plane_l = edge_l->getSupportingPlane();
    Plane3SPtr plane_r = edge_r->getSupportingPlane();
    Plane3SPtr plane_bisector;
    if (isReflex(vertex)) {
        plane_bisector = KernelWrapper::bisector(
                plane_l, KernelWrapper::opposite(plane_r));
    } else {
        plane_bisector = KernelWrapper::bisector(
                KernelWrapper::opposite(plane_l), plane_r);
    }
    Sphere3SPtr sphere = vertex->getPolygon()->getSphere();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    Vector3SPtr dir_node = KernelFactory::createVector3(
            *(node->getPoint()) - *p_center);
    Vector3SPtr normal = KernelFactory::createVector3(plane_bisector);
    Vector3SPtr dir_arc = KernelWrapper::normalize(
            KernelWrapper::cross(normal, dir_node));
    result = CircularArc::create(node, dir_arc);
    result->setEdgeLeft(getEdgeOrigin(edge_l));
    result->setEdgeRight(getEdgeOrigin(edge_r));
    result->setSupportingPlane(plane_bisector);
    data->setArc(result);
    return result;
}


bool TransSimpleSphericalSkel::init(SphericalPolygonSPtr polygon) {
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
            skel_result_->addNode(node);
            skel_result_->addArc(arc);
        } else {
            result = false;
        }
    }
    return result;
}


Point3SPtr TransSimpleSphericalSkel::offsetPoint(CircularVertexSPtr vertex, double offset) {
    Point3SPtr result = Point3SPtr();
    Sphere3SPtr sphere = vertex->getPolygon()->getSphere();
    Plane3SPtr plane_in = vertex->getEdgeIn()->getSupportingPlane();
    Plane3SPtr plane_out = vertex->getEdgeOut()->getSupportingPlane();
    Plane3SPtr plane_in_offset = KernelWrapper::offsetPlane(plane_in, offset);
    Plane3SPtr plane_out_offset = KernelWrapper::offsetPlane(plane_out, offset);
    Line3SPtr line_intersect;
    if (isReflex(vertex)) {
         line_intersect = KernelWrapper::intersection(plane_out_offset, plane_in_offset);
    } else {
         line_intersect = KernelWrapper::intersection(plane_in_offset, plane_out_offset);
    }
    result = KernelWrapper::intersection(sphere, line_intersect);
    return result;
}


double TransSimpleSphericalSkel::offsetTo(CircularEdgeSPtr edge, Point3SPtr point) {
    double result = 0.0;
    Plane3SPtr plane_edge = edge->getSupportingPlane();
    result = KernelWrapper::distance(plane_edge, point);
    if (KernelWrapper::side(plane_edge, point) < 0) {
        result *= -1.0;
    }
    return result;
}


double TransSimpleSphericalSkel::angleOf(CircularVertexSPtr vertex) {
    double result = 0.0;
    CircularEdgeSPtr edge_in = vertex->getEdgeIn();
    CircularEdgeSPtr edge_out = vertex->getEdgeOut();
    Plane3SPtr plane_in = edge_in->getSupportingPlane();
    Plane3SPtr plane_out = edge_out->getSupportingPlane();
    double angle = KernelWrapper::angle(plane_in, plane_out);
    if (isReflex(vertex)) {
        result = M_PI + angle;
    } else {
        result = M_PI - angle;
    }
    return result;
}


bool TransSimpleSphericalSkel::isEdgeEvent(CircularEdgeSPtr edge, double offset) {
    bool result = false;
    Point3SPtr p_src_offset = offsetPoint(edge->getVertexSrc(), offset);
    Point3SPtr p_dst_offset = offsetPoint(edge->getVertexDst(), offset);
    if (p_src_offset && p_dst_offset) {
        // TODO: epsilon environment is not good
        double dist_edge_event = KernelWrapper::distance(p_src_offset, p_dst_offset);
        DEBUG_VAR(dist_edge_event);
        if (dist_edge_event < epsilon_) {
            result = true;
        }
    }
    return result;
}


SphericalEdgeEventSPtr TransSimpleSphericalSkel::nextEdgeEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalEdgeEventSPtr result;
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        if (!edge->getVertexSrc()->isPointValid() ||
                !edge->getVertexDst()->isPointValid()) {
            continue;
        }
        if (isTriangle(edge)) {
            // triangle event
            continue;
        }
        CircularEdgeSPtr edge_1 = edge->prev();
        CircularEdgeSPtr edge_2 = edge->next()->next();
        CircularEdgeSPtr edge_origin_1 = getEdgeOrigin(edge_1);
        CircularEdgeSPtr edge_origin_2 = getEdgeOrigin(edge_2);
        if (edge_origin_1->getSupportingPlane() == edge_origin_2->getSupportingPlane()) {
            // edge merge event
            continue;
        }
        edge_1 = edge->prev()->prev();
        edge_2 = edge->next();
        edge_origin_1 = getEdgeOrigin(edge_1);
        edge_origin_2 = getEdgeOrigin(edge_2);
        if (edge_origin_1->getSupportingPlane() == edge_origin_2->getSupportingPlane()) {
            // edge merge event
            continue;
        }
        Point3SPtr point = vanishesAt(edge);
        if (!point) {
            continue;
        }
        double offset_current = offsetTo(edge, point);
        if (!isEdgeEvent(edge, offset_current)) {
            continue;
        }
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

SphericalSplitEventSPtr TransSimpleSphericalSkel::nextSplitEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalSplitEventSPtr result;
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        if (!vertex->isPointValid()) {
            continue;
        }
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
            bool remove_edge_l = false;
            if (edge->next()->getVertexDst() == vertex) {
                remove_edge_l = true;
            }
            bool remove_edge_r = false;
            if (edge->prev()->getVertexSrc() == vertex) {
                remove_edge_r = true;
            }
            Point3SPtr point = crashAt(vertex, edge);
            if (!point) {
                continue;
            }
            double offset_current = offsetTo(edge, point);
            Point3SPtr p_offset = offsetPoint(vertex, offset_current);
            if (!p_offset) {
                continue;
            }
            double dist_split_event = KernelWrapper::distance(point, p_offset);
            DEBUG_VAR(dist_split_event);
            if (dist_split_event > epsilon_) {
                // TODO: epsilon environment is not good
                continue;
            }
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
                if (remove_edge_l) {
                    SphericalSkelVertexDataSPtr data_l =
                            std::dynamic_pointer_cast<SphericalSkelVertexData>(
                            edge->getVertexDst()->getData());
                    node->addArc(data_l->getArc());
                }
                if (remove_edge_r) {
                    SphericalSkelVertexDataSPtr data_r =
                            std::dynamic_pointer_cast<SphericalSkelVertexData>(
                            edge->getVertexSrc()->getData());
                    node->addArc(data_r->getArc());
                }
                offset_max = offset_current;
            }
        }
    }
    return result;
}

SphericalTriangleEventSPtr TransSimpleSphericalSkel::nextTriangleEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalTriangleEventSPtr result;
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        if (!edge->getVertexSrc()->isPointValid() ||
                !edge->getVertexDst()->isPointValid()) {
            continue;
        }
        if (!isTriangle(edge)) {
            // edge event
            continue;
        }
        Point3SPtr point = vanishesAt(edge);
        if (!point) {
            continue;
        }
        double offset_current = offsetTo(edge, point);
        if (!isEdgeEvent(edge, offset_current)) {
            continue;
        }
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

SphericalDblEdgeEventSPtr TransSimpleSphericalSkel::nextDblEdgeEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalDblEdgeEventSPtr result;
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    double radius = polygon->getRadius();
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge_1 = *it_e++;
        if (!edge_1->getVertexSrc()->isPointValid() ||
                !edge_1->getVertexDst()->isPointValid()) {
            continue;
        }
        CircularEdgeSPtr edge_2 = edge_1->next();
        if (edge_2->next() != edge_1) {
            continue;
        }

        CircularVertexSPtr vertex = edge_1->getVertexDst();
        double speed = 1.0/sin(angleOf(vertex)/2.0);
        Plane3SPtr plane_1 = edge_1->getSupportingPlane();
        Plane3SPtr plane_2 = edge_2->getSupportingPlane();
        Line3SPtr line = KernelWrapper::intersection(plane_1, plane_2);
        double dist = KernelWrapper::distance(line, p_center);

        double offset_current = (dist - radius) / speed;
        if (offset_max < offset_current && offset_current <= 0.0) {
            if (!result) {
                result = SphericalDblEdgeEvent::create();
            }
            result->setOffset(offset_current + offset);
            result->setEdge1(edge_1);
            result->setEdge2(edge_2);
            offset_max = offset_current;
        }
    }
    return result;
}

SphericalLeaveEventSPtr TransSimpleSphericalSkel::nextLeaveEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalLeaveEventSPtr result;
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    double radius = polygon->getRadius();
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        if (!vertex->isPointValid()) {
            continue;
        }
        if (vertex->getData()->isHighlight()) {
            // has already left
            continue;
        }

        bool is_dbl_leave_event = false;
        SphericalSkelVertexDataSPtr data_1 =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
        CircularArcSPtr arc_1 = data_1->getArc();
        std::list<CircularVertexSPtr>::iterator it_v2 = polygon->vertices().begin();
        while (it_v2 != polygon->vertices().end()) {
            CircularVertexSPtr vertex_2 = *it_v2++;
            if (vertex == vertex_2) {
                continue;
            }
            if (!vertex_2->isPointValid()) {
                continue;
            }
            if (vertex_2->getData()->isHighlight()) {
                continue;
            }
            SphericalSkelVertexDataSPtr data_2 =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_2->getData());
            CircularArcSPtr arc_2 = data_2->getArc();
            if (arc_1->getEdgeLeft()->getSupportingPlane() == arc_2->getEdgeRight()->getSupportingPlane() &&
                    arc_1->getEdgeRight()->getSupportingPlane() == arc_2->getEdgeLeft()->getSupportingPlane()) {
                is_dbl_leave_event = true;
                break;
            }
        }
        if (is_dbl_leave_event) {
            continue;
        }

        double speed = 1.0/sin(angleOf(vertex)/2.0);
        Plane3SPtr plane_in = vertex->getEdgeIn()->getSupportingPlane();
        Plane3SPtr plane_out = vertex->getEdgeOut()->getSupportingPlane();
        Line3SPtr line = KernelWrapper::intersection(plane_in, plane_out);
        double dist = KernelWrapper::distance(line, p_center);

        if (offset != 0.0) {
            Point3SPtr p_dist = KernelWrapper::projection(line, p_center);
            Plane3SPtr plane_in_offset = KernelWrapper::offsetPlane(plane_in, -1.0);
            Plane3SPtr plane_out_offset = KernelWrapper::offsetPlane(plane_out, -1.0);
            Line3SPtr line_offset = KernelWrapper::intersection(plane_in_offset, plane_out_offset);
            Point3SPtr p_offset = KernelWrapper::projection(line_offset, p_center);
            Vector3SPtr v_dir_offset = KernelFactory::createVector3(*p_offset - *p_dist);
            Vector3SPtr v_dir_center = KernelFactory::createVector3(*p_dist - *p_center);
            //if (KernelWrapper::angle(v_dir_offset, v_dir_center) > M_PI/2.0) {
            if ((*v_dir_offset * *v_dir_center) < 0.0) {
                // moves towards center
                continue;
            }
        }

        double offset_current = (dist - radius) / speed;
        if (offset_max < offset_current && offset_current <= 0.0) {
            if (!result) {
                result = SphericalLeaveEvent::create();
            }
            result->setOffset(offset_current + offset);
            result->setVertex(vertex);
            offset_max = offset_current;
        }
    }
    return result;
}

SphericalReturnEventSPtr TransSimpleSphericalSkel::nextReturnEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalReturnEventSPtr result;
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    double radius = polygon->getRadius();
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        if (vertex->isPointValid()) {
            continue;
        }
        if (vertex->getData()->isHighlight()) {
            // has already returned
            continue;
        }

        bool is_dbl_return_event = false;
        SphericalSkelVertexDataSPtr data_1 =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
        CircularArcSPtr arc_1 = data_1->getArc();
        std::list<CircularVertexSPtr>::iterator it_v2 = polygon->vertices().begin();
        while (it_v2 != polygon->vertices().end()) {
            CircularVertexSPtr vertex_2 = *it_v2++;
            if (vertex == vertex_2) {
                continue;
            }
            if (vertex_2->isPointValid()) {
                continue;
            }
            if (vertex_2->getData()->isHighlight()) {
                continue;
            }
            SphericalSkelVertexDataSPtr data_2 =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_2->getData());
            CircularArcSPtr arc_2 = data_2->getArc();
            if (arc_1->getEdgeLeft()->getSupportingPlane() == arc_2->getEdgeRight()->getSupportingPlane() &&
                    arc_1->getEdgeRight()->getSupportingPlane() == arc_2->getEdgeLeft()->getSupportingPlane()) {
                is_dbl_return_event = true;
                break;
            }
        }
        if (is_dbl_return_event) {
            continue;
        }

        double speed = 1.0/sin(angleOf(vertex)/2.0);
        Plane3SPtr plane_in = vertex->getEdgeIn()->getSupportingPlane();
        Plane3SPtr plane_out = vertex->getEdgeOut()->getSupportingPlane();
        Line3SPtr line = KernelWrapper::intersection(plane_in, plane_out);
        double dist = KernelWrapper::distance(line, p_center);

        Point3SPtr p_dist = KernelWrapper::projection(line, p_center);
        Plane3SPtr plane_in_offset = KernelWrapper::offsetPlane(plane_in, -1.0);
        Plane3SPtr plane_out_offset = KernelWrapper::offsetPlane(plane_out, -1.0);
        Line3SPtr line_offset = KernelWrapper::intersection(plane_in_offset, plane_out_offset);
        Point3SPtr p_offset = KernelWrapper::projection(line_offset, p_center);
        Vector3SPtr v_dir_offset = KernelFactory::createVector3(*p_offset - *p_dist);
        Vector3SPtr v_dir_center = KernelFactory::createVector3(*p_center - *p_dist);
        //if (KernelWrapper::angle(v_dir_offset, v_dir_center) > M_PI/2.0) {
        if ((*v_dir_offset * *v_dir_center) < 0.0) {
            // moves away from center
            continue;
        }

        double offset_current = (radius - dist) / speed;
        if (offset_max < offset_current && offset_current <= 0.0) {
            if (!result) {
                result = SphericalReturnEvent::create();
            }
            result->setOffset(offset_current + offset);
            result->setVertex(vertex);
            offset_max = offset_current;
        }
    }
    return result;
}

SphericalDblLeaveEventSPtr TransSimpleSphericalSkel::nextDblLeaveEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalDblLeaveEventSPtr result;
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    double radius = polygon->getRadius();
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularVertexSPtr>::iterator it_v1 = polygon->vertices().begin();
    while (it_v1 != polygon->vertices().end()) {
        CircularVertexSPtr vertex_1 = *it_v1++;
        if (!vertex_1->isPointValid()) {
            continue;
        }
        if (vertex_1->getData()->isHighlight()) {
            continue;
        }
        SphericalSkelVertexDataSPtr data_1 =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_1->getData());
        CircularArcSPtr arc_1 = data_1->getArc();
        std::list<CircularVertexSPtr>::iterator it_v2 = it_v1;
        while (it_v2 != polygon->vertices().end()) {
            CircularVertexSPtr vertex_2 = *it_v2++;
            if (!vertex_2->isPointValid()) {
                continue;
            }
            if (vertex_1->next() == vertex_2 && vertex_2->next() == vertex_1) {
                // double edge event
                continue;
            }
            if (vertex_2->getData()->isHighlight()) {
                continue;
            }
            SphericalSkelVertexDataSPtr data_2 =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_2->getData());
            CircularArcSPtr arc_2 = data_2->getArc();
            if (arc_1->getEdgeLeft()->getSupportingPlane() == arc_2->getEdgeRight()->getSupportingPlane() &&
                    arc_1->getEdgeRight()->getSupportingPlane() == arc_2->getEdgeLeft()->getSupportingPlane()) {
                double speed = 1.0/sin(angleOf(vertex_1)/2.0);
                Plane3SPtr plane_in = vertex_1->getEdgeIn()->getSupportingPlane();
                Plane3SPtr plane_out = vertex_1->getEdgeOut()->getSupportingPlane();
                Line3SPtr line = KernelWrapper::intersection(plane_in, plane_out);
                double dist = KernelWrapper::distance(line, p_center);

                if (offset != 0.0) {
                    Point3SPtr p_dist = KernelWrapper::projection(line, p_center);
                    Plane3SPtr plane_in_offset = KernelWrapper::offsetPlane(plane_in, -1.0);
                    Plane3SPtr plane_out_offset = KernelWrapper::offsetPlane(plane_out, -1.0);
                    Line3SPtr line_offset = KernelWrapper::intersection(plane_in_offset, plane_out_offset);
                    Point3SPtr p_offset = KernelWrapper::projection(line_offset, p_center);
                    Vector3SPtr v_dir_offset = KernelFactory::createVector3(*p_offset - *p_dist);
                    Vector3SPtr v_dir_center = KernelFactory::createVector3(*p_dist - *p_center);
                    //if (KernelWrapper::angle(v_dir_offset, v_dir_center) > M_PI/2.0) {
                    if ((*v_dir_offset * *v_dir_center) < 0.0) {
                        // moves towards center
                        continue;
                    }
                }

                double offset_current = (dist - radius) / speed;
                if (offset_max < offset_current && offset_current <= 0.0) {
                    if (!result) {
                        result = SphericalDblLeaveEvent::create();
                    }
                    result->setOffset(offset_current + offset);
                    result->setVertex1(vertex_1);
                    result->setVertex2(vertex_2);
                    offset_max = offset_current;
                }
            }
        }
    }
    return result;
}

SphericalDblReturnEventSPtr TransSimpleSphericalSkel::nextDblReturnEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalDblReturnEventSPtr result;
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    double radius = polygon->getRadius();
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularVertexSPtr>::iterator it_v1 = polygon->vertices().begin();
    while (it_v1 != polygon->vertices().end()) {
        CircularVertexSPtr vertex_1 = *it_v1++;
        if (vertex_1->isPointValid()) {
            continue;
        }
        if (vertex_1->getData()->isHighlight()) {
            continue;
        }
        SphericalSkelVertexDataSPtr data_1 =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_1->getData());
        CircularArcSPtr arc_1 = data_1->getArc();
        std::list<CircularVertexSPtr>::iterator it_v2 = it_v1;
        while (it_v2 != polygon->vertices().end()) {
            CircularVertexSPtr vertex_2 = *it_v2++;
            if (vertex_2->isPointValid()) {
                continue;
            }
            if (vertex_2->getData()->isHighlight()) {
                continue;
            }
            SphericalSkelVertexDataSPtr data_2 =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_2->getData());
            CircularArcSPtr arc_2 = data_2->getArc();
            if (arc_1->getEdgeLeft()->getSupportingPlane() == arc_2->getEdgeRight()->getSupportingPlane() &&
                    arc_1->getEdgeRight()->getSupportingPlane() == arc_2->getEdgeLeft()->getSupportingPlane()) {
                double speed = 1.0/sin(angleOf(vertex_1)/2.0);
                Plane3SPtr plane_in = vertex_1->getEdgeIn()->getSupportingPlane();
                Plane3SPtr plane_out = vertex_1->getEdgeOut()->getSupportingPlane();
                Line3SPtr line = KernelWrapper::intersection(plane_in, plane_out);
                double dist = KernelWrapper::distance(line, p_center);

                Point3SPtr p_dist = KernelWrapper::projection(line, p_center);
                Plane3SPtr plane_in_offset = KernelWrapper::offsetPlane(plane_in, -1.0);
                Plane3SPtr plane_out_offset = KernelWrapper::offsetPlane(plane_out, -1.0);
                Line3SPtr line_offset = KernelWrapper::intersection(plane_in_offset, plane_out_offset);
                Point3SPtr p_offset = KernelWrapper::projection(line_offset, p_center);
                Vector3SPtr v_dir_offset = KernelFactory::createVector3(*p_offset - *p_dist);
                Vector3SPtr v_dir_center = KernelFactory::createVector3(*p_center - *p_dist);
                //if (KernelWrapper::angle(v_dir_offset, v_dir_center) > M_PI/2.0) {
                    if ((*v_dir_offset * *v_dir_center) < 0.0) {
                    // moves away from center
                    continue;
                }

                double offset_current = (radius - dist) / speed;
                if (offset_max < offset_current && offset_current <= 0.0) {
                    if (!result) {
                        result = SphericalDblReturnEvent::create();
                    }
                    result->setOffset(offset_current + offset);
                    result->setVertex1(vertex_1);
                    result->setVertex2(vertex_2);
                    offset_max = offset_current;
                }
            }
        }
    }
    return result;
}

SphericalVertexEventSPtr TransSimpleSphericalSkel::nextVertexEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalVertexEventSPtr result;
    Sphere3SPtr sphere = polygon->getSphere();
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularVertexSPtr>::iterator it_v1 = polygon->vertices().begin();
    while (it_v1 != polygon->vertices().end()) {
        CircularVertexSPtr vertex_1 = *it_v1++;
        if (!isReflex(vertex_1)) {
            continue;
        }
        CircularEdgeSPtr edge_1_in = vertex_1->getEdgeIn();
        if (isTriangle(edge_1_in)) {
            // triangle event
            continue;
        }
        CircularEdgeSPtr edge_1_out = vertex_1->getEdgeOut();
        Plane3SPtr plane_1_in_origin = getEdgeOrigin(edge_1_in)->getSupportingPlane();
        Plane3SPtr plane_1_out_origin = getEdgeOrigin(edge_1_out)->getSupportingPlane();
        SphericalSkelVertexDataSPtr data_1 =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_1->getData());
        CircularArcSPtr arc_1 = data_1->getArc();
        std::list<CircularVertexSPtr>::iterator it_v2 = it_v1;
        while (it_v2 != polygon->vertices().end()) {
            CircularVertexSPtr vertex_2 = *it_v2++;
            if (!isReflex(vertex_2)) {
                continue;
            }
            CircularEdgeSPtr edge_2_in = vertex_2->getEdgeIn();
            CircularEdgeSPtr edge_2_out = vertex_2->getEdgeOut();
            Plane3SPtr plane_2_in_origin = getEdgeOrigin(edge_2_in)->getSupportingPlane();
            Plane3SPtr plane_2_out_origin = getEdgeOrigin(edge_2_out)->getSupportingPlane();
            if (!(plane_1_in_origin == plane_2_out_origin ||
                    plane_2_in_origin == plane_1_out_origin)) {
                // no vertex event
                continue;
            }
            SphericalSkelVertexDataSPtr data_2 =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_2->getData());
            CircularArcSPtr arc_2 = data_2->getArc();

            CircularEdgeSPtr edge_toremove;
            Line3SPtr line_intersect;
            if (plane_1_in_origin == plane_2_out_origin) {
                line_intersect = KernelWrapper::intersection(
                        arc_1->getSupportingPlane(), arc_2->getSupportingPlane());
                if (edge_1_out->next() == edge_2_in) {
                    edge_toremove = edge_1_out;
                }
            } else {
                line_intersect = KernelWrapper::intersection(
                        arc_2->getSupportingPlane(), arc_1->getSupportingPlane());
                if (edge_2_out->next() == edge_1_in) {
                    edge_toremove = edge_2_out;
                }
            }
            if (!line_intersect) {
                continue;
            }
            Point3SPtr point = KernelWrapper::intersection(sphere, line_intersect);
            if (!point) {
                continue;
            }
            double offset_current = offsetTo(edge_1_in, point);
            if (offset_max < offset_current && offset_current <= 0.0) {
                CircularNodeSPtr node;
                if (!result) {
                    node = CircularNode::create(point);
                    result = SphericalVertexEvent::create();
                    result->setNode(node);
                }
                node = result->getNode();
                node->clear();
                node->setOffset(offset_current + offset);
                node->setPoint(point);
                result->setVertex1(vertex_1);
                result->setVertex2(vertex_2);
                node->addArc(arc_1);
                node->addArc(arc_2);
                offset_max = offset_current;
            }
        }
    }
    return result;
}

SphericalEdgeMergeEventSPtr TransSimpleSphericalSkel::nextEdgeMergeEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalEdgeMergeEventSPtr result;
    double offset_max = -std::numeric_limits<double>::max();
    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        if (!edge->getVertexSrc()->isPointValid() ||
                !edge->getVertexDst()->isPointValid()) {
            continue;
        }
        if (isTriangle(edge)) {
            // triangle event
            continue;
        }

        CircularEdgeSPtr edge_1 = edge->prev();
        CircularEdgeSPtr edge_2 = edge->next()->next();
        CircularEdgeSPtr edge_origin_1 = getEdgeOrigin(edge_1);
        CircularEdgeSPtr edge_origin_2 = getEdgeOrigin(edge_2);
        if (edge_origin_1->getSupportingPlane() != edge_origin_2->getSupportingPlane()) {
            continue;
        }

        Point3SPtr point = vanishesAt(edge);
        if (!point) {
            continue;
        }
        double offset_current = offsetTo(edge, point);
        if (!isEdgeEvent(edge, offset_current)) {
            continue;
        }
        if (offset_max < offset_current && offset_current <= 0.0) {
            CircularNodeSPtr node;
            if (!result) {
                node = CircularNode::create(point);
                result = SphericalEdgeMergeEvent::create();
                result->setNode(node);
            }
            node = result->getNode();
            node->clear();
            node->setOffset(offset_current + offset);
            node->setPoint(point);
            result->setEdge1(edge_1);
            result->setEdge2(edge_2);
            edge = edge_1;
            for (unsigned int i = 0; i < 3; i++) {
                SphericalSkelVertexDataSPtr data =
                    std::dynamic_pointer_cast<SphericalSkelVertexData>(
                    edge->getVertexDst()->getData());
                node->addArc(data->getArc());
                edge = edge->next();
            }
            offset_max = offset_current;
        }
    }
    return result;
}

SphericalInversionEventSPtr TransSimpleSphericalSkel::nextInversionEvent(SphericalPolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    SphericalInversionEventSPtr result = SphericalInversionEvent::create();
    double radius = polygon->getRadius();
    int times = (int)ceil(offset / radius) - 1;
    if ((times % 2) == 0) {
        times -= 1;
    }
    result->setOffset(times * radius);
    result->setPolygon(polygon);
    return result;
}


SphericalAbstractEventSPtr TransSimpleSphericalSkel::nextEvent(SphericalPolygonSPtr polygon, double offset) {
    SphericalAbstractEventSPtr result = SphericalAbstractEventSPtr();
    if (!polygon) {
        return result;
    }
    if (polygon->edges().size() == 0) {
        return result;
    }
    SphericalAbstractEventSPtr events[12];
    for (unsigned int i = 0; i < 12; i++) {
        events[i] = SphericalAbstractEventSPtr();
    }
    double const_offset = util::Configuration::getInstance()->getDouble(
            "algo_3d_TransSimpleSphericalSkel", "const_offset");
    if (const_offset != 0.0) {
        events[0] = SphericalConstOffsetEvent::create(offset + const_offset);
    }
    events[1] = nextEdgeEvent(polygon, offset);
    events[2] = nextSplitEvent(polygon, offset);
    events[3] = nextTriangleEvent(polygon, offset);
    events[4] = nextDblEdgeEvent(polygon, offset);
    events[5] = nextLeaveEvent(polygon, offset);
    events[6] = nextReturnEvent(polygon, offset);
    events[7] = nextDblLeaveEvent(polygon, offset);
    events[8] = nextDblReturnEvent(polygon, offset);
    events[9] = nextVertexEvent(polygon, offset);
    events[10] = nextEdgeMergeEvent(polygon, offset);
    events[11] = nextInversionEvent(polygon, offset);
    for (unsigned int i = 0; i < 12; i++) {
        if (events[i]) {
            DEBUG_VAR(events[i]->toString());
            if (!result) {
                result = events[i];
            } else if (result->getOffset() < events[i]->getOffset()) {
                result = events[i];
            } else if (result->getOffset() == events[i]->getOffset()) {
                DEBUG_VAL("WARNING: More than one possible next event.");
                DEBUG_VAL(result->toString());
                DEBUG_VAL(events[i]->toString());
                result = events[i];
            }
        }
    }
    result->setHighlight(true);
    return result;
}


SphericalPolygonSPtr TransSimpleSphericalSkel::shiftEdges(SphericalPolygonSPtr polygon, double offset) {
    Sphere3SPtr sphere = polygon->getSphere();
    SphericalPolygonSPtr result = SphericalPolygon::create(sphere);

    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        CircularVertexSPtr vertex_offset;
        Point3SPtr p_offset = offsetPoint(vertex, offset);
        if (p_offset) {
            vertex_offset = CircularVertex::create(p_offset);
        } else {
            DEBUG_VAR(vertex->toString());
            vertex_offset = CircularVertex::create(vertex->getPoint());
            vertex_offset->setPointValid(false);
        }
        SphericalSkelVertexDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
        SphericalSkelVertexDataSPtr data_offset =
                SphericalSkelVertexData::create(vertex_offset);
        data_offset->setArc(data->getArc());
        data_offset->setReflex(data->isReflex());
        data_offset->setInvert(data->doInvert());
        data->setOffsetVertex(vertex_offset);
        result->addVertex(vertex_offset);
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
        CircularVertexSPtr vertex_src_offset = data_src->getOffsetVertex();
        CircularVertexSPtr vertex_dst_offset = data_dst->getOffsetVertex();
        if (vertex_src_offset && vertex_dst_offset) {
            CircularEdgeSPtr edge_offset = CircularEdge::create(
                    vertex_src_offset, vertex_dst_offset);
            Plane3SPtr plane_edge = edge->getSupportingPlane();
            Plane3SPtr plane_edge_offset =
                    KernelWrapper::offsetPlane(plane_edge, offset);
            edge_offset->setSupportingPlane(plane_edge_offset);
            data->setOffsetEdge(edge_offset);
            result->addEdge(edge_offset);
        } else {
            DEBUG_VAR(edge->toString());
        }
    }

    return result;
}


void TransSimpleSphericalSkel::handleEdgeEvent(SphericalEdgeEventSPtr event, SphericalPolygonSPtr polygon) {
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
    vertex_data->setReflex(false);  // TODO: check
    vertex_data->setNode(node);
    CircularArcSPtr arc = createArc(vertex);

    // save and return result
    skel_result_->addArc(arc);
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void TransSimpleSphericalSkel::handleSplitEvent(SphericalSplitEventSPtr event, SphericalPolygonSPtr polygon) {
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

    bool remove_edge_l = false;
    if (edge->next()->getVertexDst() == vertex) {
        remove_edge_l = true;
    }
    bool remove_edge_r = false;
    if (edge->prev()->getVertexSrc() == vertex) {
        remove_edge_r = true;
    }

    CircularVertexSPtr vertex_l;
    CircularVertexSPtr vertex_r;
    SphericalSkelVertexDataSPtr vertex_data_l;
    SphericalSkelVertexDataSPtr vertex_data_r;
    if (!remove_edge_l && !remove_edge_r) {
        CircularVertexSPtr vertex_dst = edge->getVertexDst();
        CircularEdgeSPtr edge_l = vertex->getEdgeIn();
        edge->setVertexDst(vertex);
        vertex->setEdgeIn(edge);
        vertex_r = vertex;
        vertex_data_r = std::dynamic_pointer_cast<SphericalSkelVertexData>(
                vertex_r->getData());
        vertex_l = CircularVertex::create(node->getPoint());
        vertex_data_l = SphericalSkelVertexData::create(vertex_l);
        polygon->addVertex(vertex_l);
        edge_l->setVertexDst(vertex_l);
        vertex_l->setEdgeIn(edge_l);
        CircularEdgeSPtr edge_2 = CircularEdge::create(vertex_l, vertex_dst);
        edge_2->setSupportingPlane(edge->getSupportingPlane());
        SphericalSkelEdgeDataSPtr edge_data_2 = SphericalSkelEdgeData::create(edge_2);
        edge_data_2->setSpeed(edge_data->getSpeed());
        edge_data_2->setRotationAxis(edge_data->getRotationAxis());
        edge_data_2->setFacetOrigin(edge_data->getFacetOrigin());
        polygon->addEdge(edge_2);
    } else if (remove_edge_l && !remove_edge_r) {
        CircularEdgeSPtr edge_l = vertex->getEdgeIn();
        CircularVertexSPtr vertex_remove = edge_l->getVertexSrc();
        polygon->removeEdge(edge_l);
        vertex_remove->setEdgeIn(CircularEdgeSPtr());
        polygon->removeVertex(vertex_remove);
        edge->setVertexDst(vertex);
        vertex->setEdgeIn(edge);
        vertex_r = vertex;
        vertex_data_r = std::dynamic_pointer_cast<SphericalSkelVertexData>(
                vertex_r->getData());
    } else if (!remove_edge_l && remove_edge_r) {
        CircularEdgeSPtr edge_r = vertex->getEdgeOut();
        CircularVertexSPtr vertex_remove = edge_r->getVertexDst();
        polygon->removeEdge(edge_r);
        vertex_remove->setEdgeOut(CircularEdgeSPtr());
        polygon->removeVertex(vertex_remove);
        edge->setVertexSrc(vertex);
        vertex->setEdgeOut(edge);
        vertex_l = vertex;
        vertex_data_l = std::dynamic_pointer_cast<SphericalSkelVertexData>(
                vertex_l->getData());
    } else if (remove_edge_l && remove_edge_r) {
        CircularEdgeSPtr edges_toremove[3];
        CircularVertexSPtr vertices_toremove[3];
        edges_toremove[0] = edge;
        for (unsigned int i = 0; i < 2; i++) {
            edges_toremove[i+1] = edges_toremove[i]->next();
        }
        vertices_toremove[0] = vertex;
        for (unsigned int i = 0; i < 2; i++) {
            vertices_toremove[i+1] = vertices_toremove[i]->next();
        }
        for (unsigned int i = 0; i < 3; i++) {
            polygon->removeEdge(edges_toremove[i]);
        }
        for (unsigned int i = 0; i < 3; i++) {
            polygon->removeVertex(vertices_toremove[i]);
        }
    }

    // create arcs
    if (!remove_edge_l) {
        vertex_data_l->setReflex(false);
        vertex_data_l->setNode(node);
        CircularArcSPtr arc_l = createArc(vertex_l);
        skel_result_->addArc(arc_l);
        if (!vertex_l->getEdgeOut()->getVertexDst()->isPointValid()) {
            vertex_data_l->setInvert(false);
        }
    }
    if (!remove_edge_r) {
        vertex_data_r->setReflex(false);
        vertex_data_r->setNode(node);
        CircularArcSPtr arc_r = createArc(vertex_r);
        skel_result_->addArc(arc_r);
        if (!vertex_r->getEdgeIn()->getVertexSrc()->isPointValid()) {
            vertex_data_r->setInvert(false);
        }
    }

    // save and return result
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void TransSimpleSphericalSkel::handleTriangleEvent(SphericalTriangleEventSPtr event, SphericalPolygonSPtr polygon) {
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

void TransSimpleSphericalSkel::handleDblEdgeEvent(SphericalDblEdgeEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    SphericalSkelEdgeDataSPtr edge_data_1 =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(
            event->getEdge1()->getData());
    CircularEdgeSPtr edge_1 = edge_data_1->getOffsetEdge();
    SphericalSkelEdgeDataSPtr edge_data_2 =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(
            event->getEdge2()->getData());
    CircularEdgeSPtr edge_2 = edge_data_2->getOffsetEdge();
    CircularVertexSPtr vertex_1 = edge_1->getVertexSrc();
    CircularVertexSPtr vertex_2 = edge_1->getVertexDst();
    SphericalSkelVertexDataSPtr vertex_data_1_offset =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(
            vertex_1->getData());
    SphericalSkelVertexDataSPtr vertex_data_2_offset =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(
            vertex_2->getData());
    CircularArcSPtr arc_1 = vertex_data_1_offset->getArc();
    CircularArcSPtr arc_2 = vertex_data_2_offset->getArc();
    CircularNodeSPtr node_1 = arc_1->getNodeSrc();
    CircularNodeSPtr node_2 = arc_2->getNodeSrc();
    if (node_1->getOffset() >= node_2->getOffset()) {
        node_2->removeArc(arc_2);
        skel_result_->removeArc(arc_2);
        arc_1->setNodeDst(node_2);
        node_2->addArc(arc_1);
    } else {
        node_1->removeArc(arc_1);
        skel_result_->removeArc(arc_1);
        arc_2->setNodeDst(node_1);
        node_1->addArc(arc_2);
    }

    polygon->removeEdge(edge_1);
    polygon->removeEdge(edge_2);
    polygon->removeVertex(vertex_1);
    polygon->removeVertex(vertex_2);

    // save and return result
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void TransSimpleSphericalSkel::handleLeaveEvent(SphericalLeaveEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularVertexSPtr vertex = event->getVertex();
    SphericalSkelVertexDataSPtr data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
    CircularVertexSPtr vertex_offset = data->getOffsetVertex();
    SphericalSkelVertexDataSPtr data_offset =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_offset->getData());
    data->setHighlight(false);
    data_offset->setHighlight(true);

    Plane3SPtr plane_in = vertex_offset->getEdgeIn()->getSupportingPlane();
    Plane3SPtr plane_out = vertex_offset->getEdgeOut()->getSupportingPlane();
    Line3SPtr line = KernelWrapper::intersection(plane_in, plane_out);
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    Point3SPtr p_result = KernelWrapper::projection(line, p_center);
    vertex_offset->setPoint(p_result);
    vertex_offset->setPointValid(true);

    // save and return result
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void TransSimpleSphericalSkel::handleReturnEvent(SphericalReturnEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularVertexSPtr vertex = event->getVertex();
    SphericalSkelVertexDataSPtr data =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
    CircularVertexSPtr vertex_offset = data->getOffsetVertex();
    SphericalSkelVertexDataSPtr data_offset =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_offset->getData());
    data->setHighlight(false);
    data_offset->setHighlight(true);

    Plane3SPtr plane_in = vertex_offset->getEdgeIn()->getSupportingPlane();
    Plane3SPtr plane_out = vertex_offset->getEdgeOut()->getSupportingPlane();
    Line3SPtr line = KernelWrapper::intersection(plane_in, plane_out);
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    Point3SPtr p_result = KernelWrapper::projection(line, p_center);
    vertex_offset->setPoint(p_result);
    vertex_offset->setPointValid(true);

    // save and return result
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void TransSimpleSphericalSkel::handleDblLeaveEvent(SphericalDblLeaveEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularVertexSPtr vertex_1 = event->getVertex1();
    SphericalSkelVertexDataSPtr data_1 =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_1->getData());
    CircularVertexSPtr vertex_1_offset = data_1->getOffsetVertex();
    SphericalSkelVertexDataSPtr data_1_offset =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_1_offset->getData());
    data_1_offset->setHighlight(true);
    CircularVertexSPtr vertex_2 = event->getVertex2();
    SphericalSkelVertexDataSPtr data_2 =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_2->getData());
    CircularVertexSPtr vertex_2_offset = data_2->getOffsetVertex();
    SphericalSkelVertexDataSPtr data_2_offset =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_2_offset->getData());
    data_2_offset->setHighlight(true);

    Plane3SPtr plane_in = vertex_1_offset->getEdgeIn()->getSupportingPlane();
    Plane3SPtr plane_out = vertex_1_offset->getEdgeOut()->getSupportingPlane();
    Line3SPtr line = KernelWrapper::intersection(plane_in, plane_out);
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    Point3SPtr p_result = KernelWrapper::projection(line, p_center);
    vertex_1_offset->setPoint(p_result);
    vertex_1_offset->setPointValid(true);
    vertex_2_offset->setPoint(p_result);
    vertex_2_offset->setPointValid(true);

    // flip vertices
    //data_1_offset->setReflex(false);
    //data_2_offset->setReflex(false);

    // the following code handles a double leave event as a crash of vertices
    if (vertex_1_offset->next() == vertex_2_offset) {
        polygon->removeEdge(vertex_1_offset->getEdgeOut());
    } else {
        CircularVertexSPtr vertex_src = vertex_2_offset->getEdgeIn()->getVertexSrc();
        polygon->removeEdge(vertex_2_offset->getEdgeIn());
        vertex_1_offset->getEdgeOut()->setVertexSrc(vertex_src);
        vertex_src->setEdgeOut(vertex_1_offset->getEdgeOut());
        vertex_1_offset->setEdgeOut(CircularEdgeSPtr());
    }
    if (vertex_2_offset->next() == vertex_1_offset) {
        polygon->removeEdge(vertex_2_offset->getEdgeOut());
    } else {
        CircularVertexSPtr vertex_src = vertex_1_offset->getEdgeIn()->getVertexSrc();
        polygon->removeEdge(vertex_1_offset->getEdgeIn());
        vertex_2_offset->getEdgeOut()->setVertexSrc(vertex_src);
        vertex_src->setEdgeOut(vertex_2_offset->getEdgeOut());
        vertex_2_offset->setEdgeOut(CircularEdgeSPtr());
    }
    polygon->removeVertex(vertex_1_offset);
    polygon->removeVertex(vertex_2_offset);

    CircularArcSPtr arc_1 = data_1_offset->getArc();
    CircularArcSPtr arc_2 = data_2_offset->getArc();
    CircularNodeSPtr node_dst = arc_2->getNodeSrc();
    arc_1->setNodeDst(node_dst);
    node_dst->addArc(arc_1);
    node_dst->removeArc(arc_2);
    skel_result_->removeArc(arc_2);
    // end of crash code

    // save and return result
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void TransSimpleSphericalSkel::handleDblReturnEvent(SphericalDblReturnEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularVertexSPtr vertex_1 = event->getVertex1();
    SphericalSkelVertexDataSPtr data_1 =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_1->getData());
    CircularVertexSPtr vertex_1_offset = data_1->getOffsetVertex();
    SphericalSkelVertexDataSPtr data_1_offset =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_1_offset->getData());
    data_1_offset->setHighlight(true);
    CircularVertexSPtr vertex_2 = event->getVertex2();
    SphericalSkelVertexDataSPtr data_2 =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_2->getData());
    CircularVertexSPtr vertex_2_offset = data_2->getOffsetVertex();
    SphericalSkelVertexDataSPtr data_2_offset =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_2_offset->getData());
    data_2_offset->setHighlight(true);

    Plane3SPtr plane_in = vertex_1_offset->getEdgeIn()->getSupportingPlane();
    Plane3SPtr plane_out = vertex_1_offset->getEdgeOut()->getSupportingPlane();
    Line3SPtr line = KernelWrapper::intersection(plane_in, plane_out);
    Point3SPtr p_center = KernelFactory::createPoint3(polygon->getSphere());
    Point3SPtr p_result = KernelWrapper::projection(line, p_center);
    vertex_1_offset->setPoint(p_result);
    vertex_1_offset->setPointValid(true);
    vertex_2_offset->setPoint(p_result);
    vertex_2_offset->setPointValid(true);

    // save and return result
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void TransSimpleSphericalSkel::handleVertexEvent(SphericalVertexEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularNodeSPtr node = event->getNode();
    appendEventNode(node);

    SphericalSkelVertexDataSPtr data_1 =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(event->getVertex1()->getData());
    CircularVertexSPtr vertex_1 = data_1->getOffsetVertex();
    SphericalSkelVertexDataSPtr data_2 =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(event->getVertex2()->getData());
    CircularVertexSPtr vertex_2 = data_2->getOffsetVertex();

    CircularEdgeSPtr edge_1_in = vertex_1->getEdgeIn();
    CircularEdgeSPtr edge_1_out = vertex_1->getEdgeOut();
    Plane3SPtr plane_1_in_origin = getEdgeOrigin(edge_1_in)->getSupportingPlane();
    Plane3SPtr plane_1_out_origin = getEdgeOrigin(edge_1_out)->getSupportingPlane();
    CircularEdgeSPtr edge_2_in = vertex_2->getEdgeIn();
    CircularEdgeSPtr edge_2_out = vertex_2->getEdgeOut();
    Plane3SPtr plane_2_in_origin = getEdgeOrigin(edge_2_in)->getSupportingPlane();
    Plane3SPtr plane_2_out_origin = getEdgeOrigin(edge_2_out)->getSupportingPlane();

    vertex_1->setPoint(node->getPoint());
    vertex_2->setPoint(node->getPoint());
    CircularEdgeSPtr edge_tomerge_1;
    CircularEdgeSPtr edge_tomerge_2;
    if (plane_1_in_origin == plane_2_out_origin) {
        edge_tomerge_1 = edge_1_in;
        edge_tomerge_2 = edge_2_out;
    } else if (plane_2_in_origin == plane_1_out_origin) {
        edge_tomerge_1 = edge_2_in;
        edge_tomerge_2 = edge_1_out;
    }

    vertex_1 = edge_tomerge_1->getVertexDst();
    vertex_2 = edge_tomerge_2->getVertexSrc();
    CircularEdgeSPtr edge_1 = vertex_1->getEdgeOut();
    CircularEdgeSPtr edge_2 = vertex_2->getEdgeIn();
    edge_tomerge_1->setVertexDst(edge_tomerge_2->getVertexDst());
    edge_tomerge_1->getVertexDst()->setEdgeIn(edge_tomerge_1);
    polygon->removeEdge(edge_tomerge_2);
    edge_2->setVertexDst(vertex_1);
    vertex_1->setEdgeIn(edge_2);
    vertex_2->setEdgeIn(CircularEdgeSPtr());
    vertex_2->setEdgeOut(CircularEdgeSPtr());
    polygon->removeVertex(vertex_2);

    SphericalSkelVertexDataSPtr data_1_offset =
            std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex_1->getData());
    data_1_offset->setReflex(false);  // TODO: check
    data_1_offset->setNode(node);

    CircularArcSPtr arc = createArc(vertex_1);
    skel_result_->addArc(arc);

    // save and return result
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void TransSimpleSphericalSkel::handleEdgeMergeEvent(SphericalEdgeMergeEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    CircularNodeSPtr node = event->getNode();
    appendEventNode(node);

    SphericalSkelEdgeDataSPtr data_1 =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(event->getEdge1()->getData());
    CircularEdgeSPtr edge_1 = data_1->getOffsetEdge();
    SphericalSkelEdgeDataSPtr data_2 =
            std::dynamic_pointer_cast<SphericalSkelEdgeData>(event->getEdge2()->getData());
    CircularEdgeSPtr edge_2 = data_2->getOffsetEdge();
    CircularEdgeSPtr edge_toremove_1 = edge_1->next();
    CircularEdgeSPtr edge_toremove_2 = edge_2->prev();
    CircularVertexSPtr vertex_toremove = edge_toremove_1->getVertexDst();
    CircularVertexSPtr vertex_toremove_1 = edge_1->getVertexDst();
    CircularVertexSPtr vertex_toremove_2 = edge_2->getVertexSrc();
    CircularVertexSPtr vertex_dst = edge_2->getVertexDst();
    polygon->removeEdge(edge_toremove_1);
    polygon->removeEdge(edge_toremove_2);
    polygon->removeVertex(vertex_toremove);
    polygon->removeEdge(edge_2);
    edge_1->getVertexDst()->setEdgeIn(CircularEdgeSPtr());
    edge_1->setVertexDst(vertex_dst);
    vertex_dst->setEdgeIn(edge_1);
    polygon->removeVertex(vertex_toremove_1);
    polygon->removeVertex(vertex_toremove_2);

    // save and return result
    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void TransSimpleSphericalSkel::handleInversionEvent(SphericalInversionEventSPtr event, SphericalPolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    std::list<CircularVertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        CircularVertexSPtr vertex = *it_v++;
        SphericalSkelVertexDataSPtr data =
                std::dynamic_pointer_cast<SphericalSkelVertexData>(vertex->getData());
        if (data->doInvert()) {
            data->setReflex(!data->isReflex());
        } else {
            CircularArcSPtr arc = data->getArc();
            arc->setSupportingPlane(KernelWrapper::opposite(
                    arc->getSupportingPlane()));
            data->setInvert(true);
        }
    }

    std::list<CircularEdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        CircularEdgeSPtr edge = *it_e++;
        Plane3SPtr plane_edge = edge->getSupportingPlane();
        Plane3SPtr plane_edge_opposite = KernelWrapper::opposite(plane_edge);
        edge->setSupportingPlane(plane_edge_opposite);
    }

    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}


} }
