/**
 * @file   data/2d/skel/SplitEvent.cpp
 * @author Gernot Walzl
 * @date   2012-02-03
 */

#include "data/2d/skel/SplitEvent.h"

#include "debug.h"
#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "data/2d/skel/Node.h"
#include "data/2d/skel/SkelVertexData.h"
#include "data/2d/skel/SkelEdgeData.h"

namespace data { namespace _2d { namespace skel {

SplitEvent::SplitEvent() {
    type_ = AbstractEvent::SPLIT_EVENT;
}

SplitEvent::~SplitEvent() {
    node_.reset();
    vertex_.reset();
    edge_.reset();
}

SplitEventSPtr SplitEvent::create() {
    SplitEventSPtr result = SplitEventSPtr(new SplitEvent());
    return result;
}

NodeSPtr SplitEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void SplitEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double SplitEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getHeight();
    }
    return result;
}

VertexSPtr SplitEvent::getVertex() const {
    DEBUG_SPTR(vertex_);
    return vertex_;
}

void SplitEvent::setVertex(VertexSPtr vertex) {
    this->vertex_ = vertex;
}

EdgeSPtr SplitEvent::getEdge() const {
    DEBUG_SPTR(edge_);
    return edge_;
}

void SplitEvent::setEdge(EdgeSPtr edge) {
    this->edge_ = edge;
}

void SplitEvent::setHighlight(bool highlight) {
    if (!vertex_->hasData()) {
        SkelVertexData::create(vertex_);
    }
    vertex_->getData()->setHighlight(highlight);
    if (!edge_->hasData()) {
        SkelEdgeData::create(edge_);
    }
    edge_->getData()->setHighlight(highlight);
}

} } }
