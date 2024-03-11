/**
 * @file   data/3d/skel/SphericalDblEdgeEvent.cpp
 * @author Gernot Walzl
 * @date   2013-02-22
 */

#include "data/3d/skel/SphericalDblEdgeEvent.h"

#include "debug.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/skel/SphericalSkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

SphericalDblEdgeEvent::SphericalDblEdgeEvent() {
    this->type_ = SphericalAbstractEvent::DBL_EDGE_EVENT;
    this->offset_ = 0.0;
}

SphericalDblEdgeEvent::~SphericalDblEdgeEvent() {
    edge_1_.reset();
    edge_2_.reset();
}

SphericalDblEdgeEventSPtr SphericalDblEdgeEvent::create() {
    SphericalDblEdgeEventSPtr result = SphericalDblEdgeEventSPtr(
            new SphericalDblEdgeEvent());
    return result;
}

double SphericalDblEdgeEvent::getOffset() const {
    return this->offset_;
}

void SphericalDblEdgeEvent::setOffset(double offset) {
    this->offset_ = offset;
}

CircularEdgeSPtr SphericalDblEdgeEvent::getEdge1() const {
    return this->edge_1_;
}

void SphericalDblEdgeEvent::setEdge1(CircularEdgeSPtr edge_1) {
    this->edge_1_ = edge_1;
}

CircularEdgeSPtr SphericalDblEdgeEvent::getEdge2() const {
    return this->edge_2_;
}

void SphericalDblEdgeEvent::setEdge2(CircularEdgeSPtr edge_2) {
    this->edge_2_ = edge_2;
}

void SphericalDblEdgeEvent::setHighlight(bool highlight) {
    if (!edge_1_->hasData()) {
        SphericalSkelEdgeData::create(edge_1_);
    }
    edge_1_->getData()->setHighlight(highlight);
    if (!edge_2_->hasData()) {
        SphericalSkelEdgeData::create(edge_2_);
    }
    edge_2_->getData()->setHighlight(highlight);
}

} } }
