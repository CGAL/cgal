/**
 * @file   data/3d/skel/SphericalDblLeaveEvent.cpp
 * @author Gernot Walzl
 * @date   2013-01-28
 */

#include "data/3d/skel/SphericalDblLeaveEvent.h"

#include "debug.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/skel/SphericalSkelVertexData.h"

namespace data { namespace _3d { namespace skel {

SphericalDblLeaveEvent::SphericalDblLeaveEvent() {
    this->type_ = SphericalAbstractEvent::DBL_LEAVE_EVENT;
    this->offset_ = 0.0;
}

SphericalDblLeaveEvent::~SphericalDblLeaveEvent() {
    vertex_1_.reset();
    vertex_2_.reset();
}

SphericalDblLeaveEventSPtr SphericalDblLeaveEvent::create() {
    SphericalDblLeaveEventSPtr result = SphericalDblLeaveEventSPtr(
            new SphericalDblLeaveEvent());
    return result;
}

double SphericalDblLeaveEvent::getOffset() const {
    return this->offset_;
}

void SphericalDblLeaveEvent::setOffset(double offset) {
    this->offset_ = offset;
}

CircularVertexSPtr SphericalDblLeaveEvent::getVertex1() const {
    return this->vertex_1_;
}

void SphericalDblLeaveEvent::setVertex1(CircularVertexSPtr vertex_1) {
    this->vertex_1_ = vertex_1;
}

CircularVertexSPtr SphericalDblLeaveEvent::getVertex2() const {
    return this->vertex_2_;
}

void SphericalDblLeaveEvent::setVertex2(CircularVertexSPtr vertex_2) {
    this->vertex_2_ = vertex_2;
}

void SphericalDblLeaveEvent::setHighlight(bool highlight) {
    if (!vertex_1_->hasData()) {
        SphericalSkelVertexData::create(vertex_1_);
    }
    vertex_1_->getData()->setHighlight(highlight);
    if (!vertex_2_->hasData()) {
        SphericalSkelVertexData::create(vertex_2_);
    }
    vertex_2_->getData()->setHighlight(highlight);
}

} } }
