/**
 * @file   data/3d/skel/SphericalReturnEvent.cpp
 * @author Gernot Walzl
 * @date   2013-01-28
 */

#include "data/3d/skel/SphericalReturnEvent.h"

#include "debug.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/skel/SphericalSkelVertexData.h"

namespace data { namespace _3d { namespace skel {

SphericalReturnEvent::SphericalReturnEvent() {
    this->type_ = SphericalAbstractEvent::RETURN_EVENT;
    this->offset_ = 0.0;
}

SphericalReturnEvent::~SphericalReturnEvent() {
    vertex_.reset();
}

SphericalReturnEventSPtr SphericalReturnEvent::create() {
    SphericalReturnEventSPtr result = SphericalReturnEventSPtr(
            new SphericalReturnEvent());
    return result;
}

double SphericalReturnEvent::getOffset() const {
    return this->offset_;
}

void SphericalReturnEvent::setOffset(double offset) {
    this->offset_ = offset;
}

CircularVertexSPtr SphericalReturnEvent::getVertex() const {
    return this->vertex_;
}

void SphericalReturnEvent::setVertex(CircularVertexSPtr vertex) {
    this->vertex_ = vertex;
}

void SphericalReturnEvent::setHighlight(bool highlight) {
    if (!vertex_->hasData()) {
        SphericalSkelVertexData::create(vertex_);
    }
    vertex_->getData()->setHighlight(highlight);
}

} } }
