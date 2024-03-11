/**
 * @file   data/3d/skel/SphericalLeaveEvent.cpp
 * @author Gernot Walzl
 * @date   2013-01-28
 */

#include "data/3d/skel/SphericalLeaveEvent.h"

#include "data/3d/CircularVertex.h"
#include "data/3d/skel/SphericalSkelVertexData.h"

namespace data { namespace _3d { namespace skel {

SphericalLeaveEvent::SphericalLeaveEvent() {
    this->type_ = SphericalAbstractEvent::LEAVE_EVENT;
    this->offset_ = 0.0;
}

SphericalLeaveEvent::~SphericalLeaveEvent() {
    vertex_.reset();
}

SphericalLeaveEventSPtr SphericalLeaveEvent::create() {
    SphericalLeaveEventSPtr result = SphericalLeaveEventSPtr(
            new SphericalLeaveEvent());
    return result;
}

double SphericalLeaveEvent::getOffset() const {
    return this->offset_;
}

void SphericalLeaveEvent::setOffset(double offset) {
    this->offset_ = offset;
}

CircularVertexSPtr SphericalLeaveEvent::getVertex() const {
    return this->vertex_;
}

void SphericalLeaveEvent::setVertex(CircularVertexSPtr vertex) {
    this->vertex_ = vertex;
}

void SphericalLeaveEvent::setHighlight(bool highlight) {
    if (!vertex_->hasData()) {
        SphericalSkelVertexData::create(vertex_);
    }
    vertex_->getData()->setHighlight(highlight);
}

} } }
