/**
 * @file   data/3d/skel/SphericalDblReturnEvent.cpp
 * @author Gernot Walzl
 * @date   2013-01-28
 */

#include "data/3d/skel/SphericalDblReturnEvent.h"

#include "debug.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/skel/SphericalSkelVertexData.h"

namespace data { namespace _3d { namespace skel {

SphericalDblReturnEvent::SphericalDblReturnEvent() {
    this->type_ = SphericalAbstractEvent::DBL_RETURN_EVENT;
    this->offset_ = 0.0;
}

SphericalDblReturnEvent::~SphericalDblReturnEvent() {
    vertex_1_.reset();
    vertex_2_.reset();
}

SphericalDblReturnEventSPtr SphericalDblReturnEvent::create() {
    SphericalDblReturnEventSPtr result = SphericalDblReturnEventSPtr(
            new SphericalDblReturnEvent());
    return result;
}

double SphericalDblReturnEvent::getOffset() const {
    return this->offset_;
}

void SphericalDblReturnEvent::setOffset(double offset) {
    this->offset_ = offset;
}

CircularVertexSPtr SphericalDblReturnEvent::getVertex1() const {
    return this->vertex_1_;
}

void SphericalDblReturnEvent::setVertex1(CircularVertexSPtr vertex_1) {
    this->vertex_1_ = vertex_1;
}

CircularVertexSPtr SphericalDblReturnEvent::getVertex2() const {
    return this->vertex_2_;
}

void SphericalDblReturnEvent::setVertex2(CircularVertexSPtr vertex_2) {
    this->vertex_2_ = vertex_2;
}

void SphericalDblReturnEvent::setHighlight(bool highlight) {
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
