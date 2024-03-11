/**
 * @file   algo/2d/FastStraightSkel.cpp
 * @author Gernot Walzl
 * @date   2015-09-28
 */

#include "algo/2d/FastStraightSkel.h"

#include "debug.h"
#include "util/Timer.h"
#include "algo/Controller.h"
#include "algo/2d/KernelWrapper.h"
#include "data/2d/Polygon.h"
#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "data/2d/skel/StraightSkeleton.h"
#include "data/2d/skel/Node.h"
#include "data/2d/skel/Arc.h"
#include "data/2d/skel/ConstOffsetEvent.h"
#include "data/2d/skel/EdgeEvent.h"
#include "data/2d/skel/SplitEvent.h"
#include "data/2d/skel/TriangleEvent.h"
#include "data/2d/skel/SkelVertexData.h"
#include "data/2d/skel/SkelEdgeData.h"
#include "util/StringFactory.h"
#include <list>

namespace algo { namespace _2d {

FastStraightSkel::FastStraightSkel(PolygonSPtr polygon) {
    polygon_ = polygon;
    skel_result_ = StraightSkeleton::create();
    skel_result_->setPolygon(polygon);
}

FastStraightSkel::FastStraightSkel(PolygonSPtr polygon, ControllerSPtr controller) {
    polygon_ = polygon;
    controller_ = controller;
    skel_result_ = StraightSkeleton::create();
    skel_result_->setPolygon(polygon);
}

FastStraightSkel::~FastStraightSkel() {
    polygon_.reset();
    controller_.reset();
    skel_result_.reset();
}

FastStraightSkelSPtr FastStraightSkel::create(PolygonSPtr polygon) {
    return FastStraightSkelSPtr(new FastStraightSkel(polygon));
}

FastStraightSkelSPtr FastStraightSkel::create(PolygonSPtr polygon, ControllerSPtr controller) {
    return FastStraightSkelSPtr(new FastStraightSkel(polygon, controller));
}

void FastStraightSkel::run() {
    if (controller_) {
        controller_->wait();
        controller_->setDispPolygon(polygon_);
        controller_->setDispSkel2d(skel_result_);
    }
    DEBUG_PRINT("== Fast Straight Skeleton 2D started ==");
    double t_start = util::Timer::now();
    PolygonSPtr polygon = polygon_;
    if (init(polygon)) {
        if (controller_) {
            controller_->wait();
        }
        double offset = 0.0;
        double offset_prev = 0.0;
        int direction = -1;
        EdgeEventSPtr event = nextEvent(polygon, offset, direction);
        while (polygon->vertices().size() > 0) {
            DEBUG_VAL("-- Next Event: " << event->toString() << " --");
            if (controller_) {
                event->setHighlight(true);
                controller_->wait();
            }
            offset = event->getOffset();
            polygon = shiftEdges(polygon, offset - offset_prev);
            handleEvent(event, polygon);
            assert(polygon->isConsistent());
            assert(skel_result_->isConsistent());
            DEBUG_PRINT("-- Finished handling Event --");
            if (controller_) {
                controller_->wait();
            }
            event = nextEvent(polygon, offset, direction);
            if (!event) {
                direction *= -1;
                DEBUG_PRINT("-- direction changed --");
                event = nextEvent(polygon, offset, direction);
            }
            offset_prev = offset;
        }
        DEBUG_PRINT("== Fast Straight Skeleton 2D finished ==");
        double time = util::Timer::now() - t_start;
        skel_result_->appendDescription("time=" +
                util::StringFactory::fromDouble(time) + "; ");
        //skel_result_->appendDescription("controller=" +
        //        util::StringFactory::fromBoolean(controller_) + "; ");
        DEBUG_VAR(skel_result_->toString());
    }
}

ThreadSPtr FastStraightSkel::startThread() {
    return ThreadSPtr(new std::thread(
            std::bind(&FastStraightSkel::run, this)));
}

NodeSPtr FastStraightSkel::createNode(VertexSPtr vertex) {
    NodeSPtr result = Node::create(vertex->getPoint());
    SkelVertexDataSPtr data;
    if (vertex->hasData()) {
        data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
    } else {
        data = SkelVertexData::create(vertex);
    }
    data->setNode(result);
    DEBUG_SPTR(result);
    return result;
}

ArcSPtr FastStraightSkel::createArc(VertexSPtr vertex) {
    ArcSPtr result;
    EdgeSPtr edge_in = vertex->getEdgeIn();
    EdgeSPtr edge_out = vertex->getEdgeOut();
    SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());

    Line2SPtr line_in = edge_in->line();
    Line2SPtr line_out = edge_out->line();
    Point2SPtr p_intersect = KernelWrapper::intersection(line_in, line_out);
    if (!p_intersect) {
        return result;
    }
    EdgeSPtr edge_left = edge_in;
    double speed_in = 1.0;
    if (edge_in->hasData()) {
        SkelEdgeDataSPtr data_in = std::dynamic_pointer_cast<SkelEdgeData>(
                edge_in->getData());
        edge_left = data_in->getEdgeOrigin();
        speed_in = data_in->getSpeed();
    }
    EdgeSPtr edge_right = edge_out;
    double speed_out = 1.0;
    if (edge_out->hasData()) {
        SkelEdgeDataSPtr data_out = std::dynamic_pointer_cast<SkelEdgeData>(
                edge_out->getData());
        edge_right = data_out->getEdgeOrigin();
        speed_out = data_out->getSpeed();
    }
    Line2SPtr line_offset_in = KernelWrapper::offsetLine(line_in, speed_in);
    Line2SPtr line_offset_out = KernelWrapper::offsetLine(line_out, speed_out);
    Point2SPtr p_offset_intersect = KernelWrapper::intersection(
            line_offset_in, line_offset_out);
    Vector2SPtr dir = KernelFactory::createVector2(*p_offset_intersect - *p_intersect);

    result = Arc::create(data->getNode(), dir);
    result->setEdgeLeft(edge_left);
    result->setEdgeRight(edge_right);
    data->setArc(result);
    DEBUG_SPTR(result);
    return result;
}

bool FastStraightSkel::init(PolygonSPtr polygon) {
    WriteLock l(polygon->mutex());
    bool result = true;
    polygon->sortEdges();
    std::list<VertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        EdgeSPtr edge_in = vertex->getEdgeIn();
        EdgeSPtr edge_out = vertex->getEdgeOut();
        if (edge_in && edge_out) {
            if (!vertex->hasData()) {
                SkelVertexData::create(vertex);
            }
            if (!edge_in->hasData()) {
                SkelEdgeData::create(edge_in);
            }
            if (!edge_out->hasData()) {
                SkelEdgeData::create(edge_out);
            }
            NodeSPtr node = createNode(vertex);
            ArcSPtr arc = createArc(vertex);
            skel_result_->addNode(node);
            skel_result_->addArc(arc);
        } else {
            result = false;
        }
    }
    return result;
}

EdgeEventSPtr FastStraightSkel::nextEvent(PolygonSPtr polygon, double offset, int direction) {
    ReadLock l(polygon->mutex());
    EdgeEventSPtr result;
    if (polygon->edges().size() < 3) {
        return result;
    }

    double height_min = std::numeric_limits<double>::max();
    std::list<EdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr vertex_src = edge->getVertexSrc();
        VertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }
        SkelVertexDataSPtr data_src = std::dynamic_pointer_cast<SkelVertexData>(
                vertex_src->getData());
        SkelVertexDataSPtr data_dst = std::dynamic_pointer_cast<SkelVertexData>(
                vertex_dst->getData());
        if (data_src && data_dst) {
            ArcSPtr arc_src = data_src->getArc();
            ArcSPtr arc_dst = data_dst->getArc();
            Point2SPtr point = KernelWrapper::intersection(
                    arc_src->line(), arc_dst->line());
            if (!point) {
                continue;
            }
            if (KernelWrapper::side(edge->line(), point) == direction) {
                continue;
            }
            SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge->getData());
            double height = KernelWrapper::distance(edge->line(), point) /
                    data->getSpeed();
            if (0.0 < height && height <= height_min) {
                if (!result) {
                    result = EdgeEvent::create();
                }
                NodeSPtr node = Node::create(point);
                if (direction < 0) {
                    node->setHeight(offset + height);
                } else if (direction > 0) {
                    node->setHeight(offset - height);
                }
                node->addArc(arc_src);
                node->addArc(arc_dst);
                result->setEdge(edge);
                result->setNode(node);
                height_min = height;
            }
        }
    }

    // check if triangle event
    if (result) {
        EdgeSPtr edge_a = result->getEdge();
        EdgeSPtr edge_b = edge_a->next();
        EdgeSPtr edge_c = edge_b->next();
        if (edge_c->next() == edge_a) {
            TriangleEventSPtr event_t = TriangleEvent::create();
            NodeSPtr node = result->getNode();
            VertexDataSPtr data = edge_c->getVertexSrc()->getData();
            ArcSPtr arc = (std::dynamic_pointer_cast<SkelVertexData>(data))->getArc();
            node->addArc(arc);
            event_t->setEdge(edge_a);
            event_t->setNode(node);
            result = event_t;
        }
    }

    return result;
}

PolygonSPtr FastStraightSkel::shiftEdges(PolygonSPtr polygon, double offset) {
    PolygonSPtr result = Polygon::create();

    std::list<VertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        EdgeSPtr edge_in = vertex->getEdgeIn();
        EdgeSPtr edge_out = vertex->getEdgeOut();
        if (edge_in && edge_out) {
            double speed_in = 1.0;
            if (edge_in->hasData()) {
                SkelEdgeDataSPtr data_in = std::dynamic_pointer_cast<SkelEdgeData>(
                        edge_in->getData());
                speed_in = data_in->getSpeed();
            }
            Line2SPtr offset_line_in =
                    KernelWrapper::offsetLine(edge_in->line(), offset*speed_in);
            double speed_out = 1.0;
            if (edge_out->hasData()) {
                SkelEdgeDataSPtr data_out = std::dynamic_pointer_cast<SkelEdgeData>(
                        edge_out->getData());
                speed_out = data_out->getSpeed();
            }
            Line2SPtr offset_line_out =
                    KernelWrapper::offsetLine(edge_out->line(), offset*speed_out);
            Point2SPtr offset_point =
                    KernelWrapper::intersection(offset_line_in, offset_line_out);
            VertexSPtr offset_vertex = Vertex::create(offset_point);
            SkelVertexDataSPtr offset_data = SkelVertexData::create(offset_vertex);
            SkelVertexDataSPtr data;
            if (vertex->hasData()) {
                data = std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
                offset_data->setArc(data->getArc());
            } else {
                data = SkelVertexData::create(vertex);
            }
            data->setOffsetVertex(offset_vertex);
            result->addVertex(offset_vertex);
        }
    }

    std::list<EdgeSPtr>::iterator it_e = polygon->edges().begin();
    while (it_e != polygon->edges().end()) {
        EdgeSPtr edge = *it_e++;
        SkelVertexDataSPtr data_src = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexSrc()->getData());
        SkelVertexDataSPtr data_dst = std::dynamic_pointer_cast<SkelVertexData>(
                edge->getVertexDst()->getData());
        if (data_src && data_dst) {
            EdgeSPtr offset_edge = Edge::create(
                    data_src->getOffsetVertex(), data_dst->getOffsetVertex());
            SkelEdgeDataSPtr offset_data = SkelEdgeData::create(offset_edge);
            SkelEdgeDataSPtr data;
            if (edge->hasData()) {
                data = std::dynamic_pointer_cast<SkelEdgeData>(edge->getData());
                offset_data->setEdgeOrigin(data->getEdgeOrigin());
                offset_data->setSpeed(data->getSpeed());
            } else {
                data = SkelEdgeData::create(edge);
            }
            data->setOffsetEdge(offset_edge);
            result->addEdge(offset_edge);
        }
    }

    return result;
}

void FastStraightSkel::appendEventNode(NodeSPtr node) {
    std::list<ArcWPtr>::iterator it_a = node->arcs().begin();
    while (it_a != node->arcs().end()) {
        ArcWPtr arc_wptr = *it_a++;
        if (!arc_wptr.expired()) {
            ArcSPtr arc = ArcSPtr(arc_wptr);
            arc->setNodeDst(node);
            arc->setNodeDstListIt(
                    std::find(node->arcs().begin(), node->arcs().end(), arc_wptr));
        }
    }
    skel_result_->addNode(node);
}

void FastStraightSkel::handleEvent(EdgeEventSPtr event, PolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    if (event->getType() == AbstractEvent::EDGE_EVENT) {
        EdgeSPtr edge = event->getEdge();
        SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                edge->getData());
        EdgeSPtr edge_result = edge_data->getOffsetEdge();
        VertexSPtr vertex_result = edge_result->getVertexSrc();

        edge_result->next()->setVertexSrc(vertex_result);
        vertex_result->setEdgeOut(edge_result->next());
        polygon->removeEdge(edge_result);
        edge_result->getVertexDst()->setEdgeOut(EdgeSPtr());
        polygon->removeVertex(edge_result->getVertexDst());

        SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                vertex_result->getData());
        data->setNode(node);
        ArcSPtr arc_out = createArc(vertex_result);
        skel_result_->addArc(arc_out);
    } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
        TriangleEventSPtr event_t = std::dynamic_pointer_cast<TriangleEvent>(event);
        EdgeSPtr edges_toremove[3];
        event_t->getEdges(edges_toremove);
        for (unsigned int i = 0; i<3; i++) {
            SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                    edges_toremove[i]->getData());
            polygon->removeEdge(edge_data->getOffsetEdge());
        }

        VertexSPtr vertices_toremove[3];
        event_t->getVertices(vertices_toremove);
        for (unsigned int i = 0; i<3; i++) {
            SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                    vertices_toremove[i]->getData());
            polygon->removeVertex(vertex_data->getOffsetVertex());
        }
    }

    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

StraightSkeletonSPtr FastStraightSkel::getResult() const {
    return this->skel_result_;
}

} }
