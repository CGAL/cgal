/**
 * @file   data/3d/skel/SphericalVertexEvent.cpp
 * @author Gernot Walzl
 * @date   2013-02-20
 */

#include "data/3d/skel/SphericalVertexEvent.h"

#include "debug.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/SphericalSkelVertexData.h"

namespace data { namespace _3d { namespace skel {

SphericalVertexEvent::SphericalVertexEvent() {
    type_ = SphericalAbstractEvent::VERTEX_EVENT;
}

SphericalVertexEvent::~SphericalVertexEvent() {
    node_.reset();
    vertex_1_.reset();
    vertex_2_.reset();
}

SphericalVertexEventSPtr SphericalVertexEvent::create() {
    SphericalVertexEventSPtr result = SphericalVertexEventSPtr(new SphericalVertexEvent());
    return result;
}

CircularNodeSPtr SphericalVertexEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void SphericalVertexEvent::setNode(CircularNodeSPtr node) {
    this->node_ = node;
}

double SphericalVertexEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

CircularVertexSPtr SphericalVertexEvent::getVertex1() const {
    DEBUG_SPTR(vertex_1_);
    return vertex_1_;
}

void SphericalVertexEvent::setVertex1(CircularVertexSPtr vertex_1) {
    this->vertex_1_ = vertex_1;
}

CircularVertexSPtr SphericalVertexEvent::getVertex2() const {
    DEBUG_SPTR(vertex_2_);
    return vertex_2_;
}

void SphericalVertexEvent::setVertex2(CircularVertexSPtr vertex_2) {
    this->vertex_2_ = vertex_2;
}

void SphericalVertexEvent::setHighlight(bool highlight) {
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
