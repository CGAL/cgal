/**
 * @file   data/3d/skel/SphericalEdgeEvent.cpp
 * @author Gernot Walzl
 * @date   2012-11-29
 */

#include "data/3d/skel/SphericalEdgeEvent.h"

#include "debug.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/SphericalSkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

SphericalEdgeEvent::SphericalEdgeEvent() {
    type_ = SphericalAbstractEvent::EDGE_EVENT;
}

SphericalEdgeEvent::~SphericalEdgeEvent() {
    node_.reset();
    edge_.reset();
}

SphericalEdgeEventSPtr SphericalEdgeEvent::create() {
    SphericalEdgeEventSPtr result = SphericalEdgeEventSPtr(new SphericalEdgeEvent());
    return result;
}

CircularNodeSPtr SphericalEdgeEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void SphericalEdgeEvent::setNode(CircularNodeSPtr node) {
    this->node_ = node;
}

double SphericalEdgeEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

CircularEdgeSPtr SphericalEdgeEvent::getEdge() const {
    DEBUG_SPTR(edge_);
    return edge_;
}

void SphericalEdgeEvent::setEdge(CircularEdgeSPtr edge) {
    this->edge_ = edge;
}

void SphericalEdgeEvent::setHighlight(bool highlight) {
    if (!edge_->hasData()) {
        SphericalSkelEdgeData::create(edge_);
    }
    edge_->getData()->setHighlight(highlight);
}

} } }
