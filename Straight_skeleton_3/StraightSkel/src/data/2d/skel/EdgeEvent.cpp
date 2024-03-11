/**
 * @file   data/2d/skel/EdgeEvent.cpp
 * @author Gernot Walzl
 * @date   2012-02-03
 */

#include "data/2d/skel/EdgeEvent.h"

#include "debug.h"
#include "data/2d/Edge.h"
#include "data/2d/skel/Node.h"
#include "data/2d/skel/SkelEdgeData.h"

namespace data { namespace _2d { namespace skel {

EdgeEvent::EdgeEvent() {
    type_ = AbstractEvent::EDGE_EVENT;
}

EdgeEvent::~EdgeEvent() {
    node_.reset();
    edge_.reset();
}

EdgeEventSPtr EdgeEvent::create() {
    EdgeEventSPtr result = EdgeEventSPtr(new EdgeEvent());
    return result;
}

NodeSPtr EdgeEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void EdgeEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double EdgeEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getHeight();
    }
    return result;
}

EdgeSPtr EdgeEvent::getEdge() const {
    DEBUG_SPTR(edge_);
    return edge_;
}

void EdgeEvent::setEdge(EdgeSPtr edge) {
    this->edge_ = edge;
}

void EdgeEvent::setHighlight(bool highlight) {
    if (!edge_->hasData()) {
        SkelEdgeData::create(edge_);
    }
    edge_->getData()->setHighlight(highlight);
}

} } }
