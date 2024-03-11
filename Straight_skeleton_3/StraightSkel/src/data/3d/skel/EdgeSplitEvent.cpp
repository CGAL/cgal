/**
 * @file   data/3d/skel/EdgeSplitEvent.cpp
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#include "data/3d/skel/EdgeSplitEvent.h"

#include "debug.h"
#include "data/3d/Edge.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

EdgeSplitEvent::EdgeSplitEvent() {
    type_ = AbstractEvent::EDGE_SPLIT_EVENT;
}

EdgeSplitEvent::~EdgeSplitEvent() {
    node_.reset();
    edge1_.reset();
    edge2_.reset();
}

EdgeSplitEventSPtr EdgeSplitEvent::create() {
    EdgeSplitEventSPtr result = EdgeSplitEventSPtr(new EdgeSplitEvent());
    return result;
}

NodeSPtr EdgeSplitEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void EdgeSplitEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double EdgeSplitEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr EdgeSplitEvent::getEdge1() const {
    DEBUG_SPTR(edge1_);
    return edge1_;
}

void EdgeSplitEvent::setEdge1(EdgeSPtr edge1) {
    this->edge1_ = edge1;
}

EdgeSPtr EdgeSplitEvent::getEdge2() const {
    DEBUG_SPTR(edge2_);
    return edge2_;
}

void EdgeSplitEvent::setEdge2(EdgeSPtr edge2) {
    this->edge2_ = edge2;
}

void EdgeSplitEvent::setHighlight(bool highlight) {
    if (!edge1_->hasData()) {
        SkelEdgeData::create(edge1_);
    }
    edge1_->getData()->setHighlight(highlight);
    if (!edge2_->hasData()) {
        SkelEdgeData::create(edge2_);
    }
    edge2_->getData()->setHighlight(highlight);
}

} } }
