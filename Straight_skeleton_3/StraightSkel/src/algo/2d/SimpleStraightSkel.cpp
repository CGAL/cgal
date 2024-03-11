/**
 * @file   algo/2d/SimpleStraightSkel.cpp
 * @author Gernot Walzl
 * @date   2012-02-06
 */

#include "algo/2d/SimpleStraightSkel.h"

#include "debug.h"
#include "util/Configuration.h"
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

namespace algo { namespace _2d {

SimpleStraightSkel::SimpleStraightSkel(PolygonSPtr polygon) {
    polygon_ = polygon;
    skel_result_ = StraightSkeleton::create();
    skel_result_->setPolygon(polygon);
}

SimpleStraightSkel::SimpleStraightSkel(PolygonSPtr polygon, ControllerSPtr controller) {
    polygon_ = polygon;
    controller_ = controller;
    skel_result_ = StraightSkeleton::create();
    skel_result_->setPolygon(polygon);
}

SimpleStraightSkel::~SimpleStraightSkel() {
    polygon_.reset();
    controller_.reset();
    skel_result_.reset();
}

SimpleStraightSkelSPtr SimpleStraightSkel::create(PolygonSPtr polygon) {
    return SimpleStraightSkelSPtr(new SimpleStraightSkel(polygon));
}

SimpleStraightSkelSPtr SimpleStraightSkel::create(PolygonSPtr polygon, ControllerSPtr controller) {
    return SimpleStraightSkelSPtr(new SimpleStraightSkel(polygon, controller));
}

void SimpleStraightSkel::run() {
    if (controller_) {
        controller_->wait();
        controller_->setDispPolygon(polygon_);
        controller_->setDispSkel2d(skel_result_);
    }
    DEBUG_PRINT("== Straight Skeleton 2D started ==");
    double t_start = util::Timer::now();
    PolygonSPtr polygon = polygon_;
    unsigned int i = 0;
    if (init(polygon)) {
        if (controller_) {
            controller_->wait();
        }
        double offset = 0.0;
        double offset_prev = 0.0;
        std::list<AbstractEventSPtr> events = nextEvent(polygon, offset);
        while (events.size() > 0) {
            std::list<AbstractEventSPtr>::iterator it_e = events.begin();
            while (it_e != events.end()) {
                AbstractEventSPtr event = *it_e++;
                DEBUG_VAL("-- Next Event: " << event->toString() << " --");
                event->setHighlight(true);
            }
            if (controller_) {
                controller_->wait();
            }
            offset = events.front()->getOffset();
            polygon = shiftEdges(polygon, offset - offset_prev);
            it_e = events.begin();
            while (it_e != events.end()) {
                AbstractEventSPtr event = *it_e++;
                if (event->getType() == AbstractEvent::CONST_OFFSET_EVENT) {
                    event->setPolygonResult(polygon);
                    skel_result_->addEvent(event);
                    bool screenshot_on_const_offset_event =
                            util::Configuration::getInstance()->getBool(
                            "algo_2d_SimpleStraightSkel", "screenshot_on_const_offset_event");
                    if (controller_ && screenshot_on_const_offset_event) {
                        controller_->screenshot();
                    }
                } else if (event->getType() == AbstractEvent::EDGE_EVENT) {
                    handleEdgeEvent(std::dynamic_pointer_cast<EdgeEvent>(event), polygon);
                } else if (event->getType() == AbstractEvent::SPLIT_EVENT) {
                    handleSplitEvent(std::dynamic_pointer_cast<SplitEvent>(event), polygon);
                } else  if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
                    handleTriangleEvent(std::dynamic_pointer_cast<TriangleEvent>(event), polygon);
                }
            }
            assert(polygon->isConsistent());
            assert(skel_result_->isConsistent());
            DEBUG_PRINT("-- Finished handling Event --");
            i++;
            DEBUG_VAR(i);
            if (controller_) {
                controller_->wait();
            }
            events = nextEvent(polygon, offset);
            offset_prev = offset;
        }
        DEBUG_PRINT("== Straight Skeleton 2D finished ==");
        double time = util::Timer::now() - t_start;
        skel_result_->appendDescription("time=" +
                util::StringFactory::fromDouble(time) + "; ");
        //skel_result_->appendDescription("controller=" +
        //        util::StringFactory::fromBoolean(controller_) + "; ");
        DEBUG_VAR(skel_result_->toString());
    }
}

ThreadSPtr SimpleStraightSkel::startThread() {
    return ThreadSPtr(new std::thread(
            std::bind(&SimpleStraightSkel::run, this)));
}


NodeSPtr SimpleStraightSkel::createNode(VertexSPtr vertex) {
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

ArcSPtr SimpleStraightSkel::createArc(VertexSPtr vertex) {
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

bool SimpleStraightSkel::init(PolygonSPtr polygon) {
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

std::list<EdgeEventSPtr> SimpleStraightSkel::nextEdgeEvent(PolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    std::list<EdgeEventSPtr> result;
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
            if (KernelWrapper::side(edge->line(), point) < 0) {
                continue;
            }
            SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(
                    edge->getData());
            double height = KernelWrapper::distance(edge->line(), point) /
                    data->getSpeed();
            if (0.0 < height && height <= height_min) {
                if (height < height_min) {
                    result.clear();
                }
                EdgeEventSPtr event = EdgeEvent::create();
                NodeSPtr node = Node::create(point);
                node->setHeight(height + offset);
                node->addArc(arc_src);
                node->addArc(arc_dst);
                event->setEdge(edge);
                event->setNode(node);
                result.push_back(event);
                height_min = height;
            }
        }
    }

    // check if triangle event
    if (result.size() > 0) {
        std::list<EdgeEventSPtr>::iterator it_e = result.begin();
        while (it_e != result.end()) {
            std::list<EdgeEventSPtr>::iterator it_current = it_e;
            EdgeEventSPtr event = *it_e++;
            if (event->getType() == AbstractEvent::EDGE_EVENT) {
                EdgeSPtr edge_a = event->getEdge();
                EdgeSPtr edge_b = edge_a->next();
                EdgeSPtr edge_c = edge_b->next();
                if (edge_c->next() == edge_a) {
                    TriangleEventSPtr event_t = TriangleEvent::create();
                    NodeSPtr node = event->getNode();
                    VertexDataSPtr data = edge_c->getVertexSrc()->getData();
                    ArcSPtr arc = (std::dynamic_pointer_cast<SkelVertexData>(data))->getArc();
                    node->addArc(arc);
                    event_t->setEdge(edge_a);
                    event_t->setNode(node);
                    result.push_back(event_t);

                    std::list<EdgeEventSPtr>::iterator it_e2 = it_e;
                    while (it_e2 != result.end()) {
                        std::list<EdgeEventSPtr>::iterator it_current_2 = it_e2;
                        EdgeEventSPtr event_2 = *it_e2++;
                        if (event_2->getEdge() == edge_b ||
                                event_2->getEdge() == edge_c) {
                            result.erase(it_current_2);
                        }
                    }
                    it_e = it_current;
                    it_e++;
                    result.erase(it_current);
                }
            }
        }
    }

    return result;
}

Point2SPtr SimpleStraightSkel::crashAt(VertexSPtr vertex, EdgeSPtr edge) {
    Point2SPtr result;
    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
            vertex->getData());
    ArcSPtr arc = vertex_data->getArc();
    Line2SPtr line_arc = arc->line();

    Point2SPtr p_vertex = vertex->getPoint();
    EdgeSPtr edge_vertex = vertex->getEdgeOut();
    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            edge_vertex->getData());
    Line2SPtr line_edge_vertex_offset = KernelWrapper::offsetLine(
            edge_vertex->line(), edge_data->getSpeed());
    Point2SPtr p_vertex_offset = KernelWrapper::intersection(
            line_arc, line_edge_vertex_offset);
    if (!p_vertex_offset) {
        return result;
    }
    double speed_vertex = KernelWrapper::distance(p_vertex_offset, p_vertex);

    Line2SPtr line_edge = edge->line();
    Point2SPtr p_edge = KernelWrapper::intersection(line_arc, line_edge);
    if (!p_edge) {  // parallel
        return result;
    }
    edge_data = std::dynamic_pointer_cast<SkelEdgeData>(edge->getData());
    Line2SPtr line_edge_offset = KernelWrapper::offsetLine(
            line_edge, edge_data->getSpeed());
    Point2SPtr p_edge_offset = KernelWrapper::intersection(
            line_arc, line_edge_offset);
    double speed_edge = KernelWrapper::distance(p_edge_offset, p_edge);

    Vector2SPtr dir = KernelFactory::createVector2(*p_edge - *p_vertex);
    //if (KernelWrapper::angle(dir, arc->getDirection()) < M_PI/2.0) {
    if ((*dir * *(arc->getDirection())) > 0.0) {
        double distance = KernelWrapper::distance(p_edge, p_vertex);
        double dist_vertex = distance * speed_vertex/(speed_vertex+speed_edge);
        Point2SPtr p_crash = KernelWrapper::offsetPoint(p_vertex, arc->getDirection(), dist_vertex);
        if (KernelWrapper::side(line_edge, p_crash) >= 0) {
            SkelVertexDataSPtr data_src = std::dynamic_pointer_cast<SkelVertexData>(
                    edge->getVertexSrc()->getData());
            SkelVertexDataSPtr data_dst = std::dynamic_pointer_cast<SkelVertexData>(
                    edge->getVertexDst()->getData());
            ArcSPtr arc_src = data_src->getArc();
            ArcSPtr arc_dst = data_dst->getArc();
            if (KernelWrapper::side(arc_src->line(), p_crash) <= 0 &&
                    KernelWrapper::side(arc_dst->line(), p_crash) >= 0) {
                result = p_crash;
            }
        }
    }
    return result;
}

std::list<SplitEventSPtr> SimpleStraightSkel::nextSplitEvent(PolygonSPtr polygon, double offset) {
    ReadLock l(polygon->mutex());
    std::list<SplitEventSPtr> result;
    double height_min = std::numeric_limits<double>::max();
    std::list<VertexSPtr>::iterator it_v = polygon->vertices().begin();
    while (it_v != polygon->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (vertex->isReflex()) {
            std::list<EdgeSPtr>::iterator it_e = polygon->edges().begin();
            while (it_e != polygon->edges().end()) {
                EdgeSPtr edge = *it_e++;
                if (edge == vertex->getEdgeIn() ||
                        edge == vertex->getEdgeOut()) {
                    continue;
                }
                if (edge == vertex->getEdgeIn()->prev() ||
                        edge == vertex->getEdgeOut()->next()) {
                    // sticking event
                    continue;
                }
                Point2SPtr p_crash = crashAt(vertex, edge);
                if (!p_crash) {
                    continue;
                }
                SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(
                        edge->getData());
                double height = KernelWrapper::distance(edge->line(), p_crash) /
                        data->getSpeed();
                if (0.0 < height && height <= height_min) {
                    if (height < height_min) {
                        result.clear();
                    }
                    NodeSPtr node = Node::create(p_crash);
                    node->setHeight(height + offset);
                    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                            vertex->getData());
                    ArcSPtr arc = vertex_data->getArc();
                    node->addArc(arc);
                    SplitEventSPtr event = SplitEvent::create();
                    event->setNode(node);
                    event->setVertex(vertex);
                    event->setEdge(edge);
                    result.push_back(event);
                    height_min = height;
                }
            }
        }
    }
    return result;
}

std::list<AbstractEventSPtr> SimpleStraightSkel::nextEvent(PolygonSPtr polygon, double offset) {
    std::list<AbstractEventSPtr> result;
    if (!polygon) {
        return result;
    }
    if (polygon->edges().size() == 0) {
        return result;
    }
    double const_offset = util::Configuration::getInstance()->getDouble(
            "algo_2d_SimpleStraightSkel", "const_offset");
    if (const_offset != 0.0) {
        double next_offset = floor(offset/const_offset + 1.0) * const_offset;
        if (next_offset <= offset) {
            next_offset += const_offset;
        }
        ConstOffsetEventSPtr event = ConstOffsetEvent::create(next_offset);
        result.push_back(event);
    }
    std::list<EdgeEventSPtr> edge_events = nextEdgeEvent(polygon, offset);
    std::list<EdgeEventSPtr>::iterator it_ee = edge_events.begin();
    while (it_ee != edge_events.end()) {
        EdgeEventSPtr event = *it_ee++;
        result.push_back(event);
    }
    std::list<SplitEventSPtr> split_events = nextSplitEvent(polygon, offset);
    std::list<SplitEventSPtr>::iterator it_se = split_events.begin();
    while (it_se != split_events.end()) {
        SplitEventSPtr event = *it_se++;
        result.push_back(event);
    }
    if (result.size() > 0) {
        double offset_min = std::numeric_limits<double>::max();
        std::list<AbstractEventSPtr>::iterator it_e = result.begin();
        while (it_e != result.end()) {
            AbstractEventSPtr event = *it_e++;
            if (event->getOffset() < offset_min) {
                offset_min = event->getOffset();
            }
        }
        it_e = result.begin();
        while (it_e != result.end()) {
            std::list<AbstractEventSPtr>::iterator it_current = it_e;
            AbstractEventSPtr event = *it_e++;
            if (event->getOffset() > offset_min) {
                result.erase(it_current);
            }
        }
    }
    return result;
}

PolygonSPtr SimpleStraightSkel::shiftEdges(PolygonSPtr polygon, double offset) {
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

void SimpleStraightSkel::appendEventNode(NodeSPtr node) {
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

void SimpleStraightSkel::handleEdgeEvent(EdgeEventSPtr event, PolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    NodeSPtr node = event->getNode();
    appendEventNode(node);

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

    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handleSplitEvent(SplitEventSPtr event, PolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
            event->getVertex()->getData());
    VertexSPtr vertex = vertex_data->getOffsetVertex();
    vertex->setPoint(node->getPoint());
    SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
            event->getEdge()->getData());
    EdgeSPtr edge = edge_data->getOffsetEdge();

    VertexSPtr vertex_left = vertex;
    VertexSPtr vertex_right = Vertex::create(vertex->getPoint());
    SkelVertexDataSPtr vertex_data_left = std::dynamic_pointer_cast<SkelVertexData>(vertex_left->getData());
    SkelVertexDataSPtr vertex_data_right = SkelVertexData::create(vertex_right);

    EdgeSPtr edge_left = edge;
    EdgeSPtr edge_right = Edge::create(edge->getVertexSrc(), vertex_right);
    SkelEdgeDataSPtr edge_data_left = std::dynamic_pointer_cast<SkelEdgeData>(edge_left->getData());
    SkelEdgeDataSPtr edge_data_right = SkelEdgeData::create(edge_right);
    edge_data_right->setEdgeOrigin(edge_data_left->getEdgeOrigin());
    edge_data_right->setSpeed(edge_data_left->getSpeed());

    vertex_right->setEdgeOut(vertex->getEdgeOut());
    vertex->getEdgeOut()->setVertexSrc(vertex_right);
    vertex_left->setEdgeOut(edge_left);
    edge_left->setVertexSrc(vertex_left);
    polygon->addVertex(vertex_right);
    polygon->addEdge(edge_right);
    polygon->sortEdges();

    vertex_data_left->setNode(node);
    vertex_data_right->setNode(node);
    ArcSPtr arc_out_left = createArc(vertex_left);
    ArcSPtr arc_out_right = createArc(vertex_right);

    event->setPolygonResult(polygon);
    skel_result_->addArc(arc_out_left);
    skel_result_->addArc(arc_out_right);
    skel_result_->addEvent(event);
}

void SimpleStraightSkel::handleTriangleEvent(TriangleEventSPtr event, PolygonSPtr polygon) {
    WriteLock l(skel_result_->mutex());

    NodeSPtr node = event->getNode();
    appendEventNode(node);

    EdgeSPtr edges_toremove[3];
    event->getEdges(edges_toremove);
    for (unsigned int i = 0; i<3; i++) {
        SkelEdgeDataSPtr edge_data = std::dynamic_pointer_cast<SkelEdgeData>(
                edges_toremove[i]->getData());
        polygon->removeEdge(edge_data->getOffsetEdge());
    }

    VertexSPtr vertices_toremove[3];
    event->getVertices(vertices_toremove);
    for (unsigned int i = 0; i<3; i++) {
        SkelVertexDataSPtr vertex_data = std::dynamic_pointer_cast<SkelVertexData>(
                vertices_toremove[i]->getData());
        polygon->removeVertex(vertex_data->getOffsetVertex());
    }

    event->setPolygonResult(polygon);
    skel_result_->addEvent(event);
}


StraightSkeletonSPtr SimpleStraightSkel::getResult() const {
    return this->skel_result_;
}

} }
