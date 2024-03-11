/**
 * @file   data/3d/skel/SurfaceEvent.cpp
 * @author Gernot Walzl
 * @date   2012-09-10
 */

#include "data/3d/skel/SurfaceEvent.h"

#include "debug.h"
#include "data/3d/Edge.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

SurfaceEvent::SurfaceEvent() {
    type_ = AbstractEvent::SURFACE_EVENT;
}

SurfaceEvent::~SurfaceEvent() {
    node_.reset();
    edge1_.reset();
    edge2_.reset();
}

SurfaceEventSPtr SurfaceEvent::create() {
    SurfaceEventSPtr result = SurfaceEventSPtr(new SurfaceEvent());
    return result;
}

NodeSPtr SurfaceEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void SurfaceEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double SurfaceEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr SurfaceEvent::getEdge1() const {
    DEBUG_SPTR(edge1_);
    return edge1_;
}

void SurfaceEvent::setEdge1(EdgeSPtr edge1) {
    this->edge1_ = edge1;
}

EdgeSPtr SurfaceEvent::getEdge2() const {
    DEBUG_SPTR(edge2_);
    return edge2_;
}

void SurfaceEvent::setEdge2(EdgeSPtr edge2) {
    this->edge2_ = edge2;
}

void SurfaceEvent::setHighlight(bool highlight) {
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
