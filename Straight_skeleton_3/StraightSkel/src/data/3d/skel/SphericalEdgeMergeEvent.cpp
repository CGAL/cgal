/**
 * @file   data/3d/skel/SphericalEdgeMergeEvent.cpp
 * @author Gernot Walzl
 * @date   2013-02-08
 */

#include "data/3d/skel/SphericalEdgeMergeEvent.h"

#include "debug.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/SphericalSkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

SphericalEdgeMergeEvent::SphericalEdgeMergeEvent() {
    type_ = SphericalAbstractEvent::EDGE_MERGE_EVENT;
}

SphericalEdgeMergeEvent::~SphericalEdgeMergeEvent() {
    node_.reset();
    edge_1_.reset();
    edge_2_.reset();
}

SphericalEdgeMergeEventSPtr SphericalEdgeMergeEvent::create() {
    SphericalEdgeMergeEventSPtr result = SphericalEdgeMergeEventSPtr(new SphericalEdgeMergeEvent());
    return result;
}

CircularNodeSPtr SphericalEdgeMergeEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void SphericalEdgeMergeEvent::setNode(CircularNodeSPtr node) {
    this->node_ = node;
}

double SphericalEdgeMergeEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

CircularEdgeSPtr SphericalEdgeMergeEvent::getEdge1() const {
    DEBUG_SPTR(edge_1_);
    return edge_1_;
}

void SphericalEdgeMergeEvent::setEdge1(CircularEdgeSPtr edge_1) {
    this->edge_1_ = edge_1;
}

CircularEdgeSPtr SphericalEdgeMergeEvent::getEdge2() const {
    DEBUG_SPTR(edge_2_);
    return edge_2_;
}

void SphericalEdgeMergeEvent::setEdge2(CircularEdgeSPtr edge_2) {
    this->edge_2_ = edge_2;
}

void SphericalEdgeMergeEvent::setHighlight(bool highlight) {
    if (!edge_1_->hasData()) {
        SphericalSkelEdgeData::create(edge_1_);
    }
    edge_1_->getData()->setHighlight(highlight);
    if (!edge_2_->hasData()) {
        SphericalSkelEdgeData::create(edge_2_);
    }
    edge_2_->getData()->setHighlight(highlight);
    CircularEdgeSPtr edge_toremove_1 = edge_1_->next();
    CircularEdgeSPtr edge_toremove_2 = edge_toremove_1->next();
    if (!edge_toremove_1->hasData()) {
        SphericalSkelEdgeData::create(edge_toremove_1);
    }
    edge_toremove_1->getData()->setHighlight(highlight);
    if (!edge_toremove_2->hasData()) {
        SphericalSkelEdgeData::create(edge_toremove_2);
    }
    edge_toremove_2->getData()->setHighlight(highlight);
}

} } }
