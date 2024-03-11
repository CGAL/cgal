/**
 * @file   data/3d/skel/PolyhedronSplitEvent.cpp
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#include "data/3d/skel/PolyhedronSplitEvent.h"

#include "debug.h"
#include "data/3d/Edge.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

PolyhedronSplitEvent::PolyhedronSplitEvent() {
    type_ = AbstractEvent::POLYHEDRON_SPLIT_EVENT;
}

PolyhedronSplitEvent::~PolyhedronSplitEvent() {
    node_.reset();
    edge1_.reset();
    edge2_.reset();
}

PolyhedronSplitEventSPtr PolyhedronSplitEvent::create() {
    PolyhedronSplitEventSPtr result = PolyhedronSplitEventSPtr(new PolyhedronSplitEvent());
    return result;
}

NodeSPtr PolyhedronSplitEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void PolyhedronSplitEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double PolyhedronSplitEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr PolyhedronSplitEvent::getEdge1() const {
    DEBUG_SPTR(edge1_);
    return edge1_;
}

void PolyhedronSplitEvent::setEdge1(EdgeSPtr edge1) {
    this->edge1_ = edge1;
}

EdgeSPtr PolyhedronSplitEvent::getEdge2() const {
    DEBUG_SPTR(edge2_);
    return edge2_;
}

void PolyhedronSplitEvent::setEdge2(EdgeSPtr edge2) {
    this->edge2_ = edge2;
}

void PolyhedronSplitEvent::setHighlight(bool highlight) {
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
