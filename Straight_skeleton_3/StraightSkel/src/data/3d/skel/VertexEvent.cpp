/**
 * @file   data/3d/skel/VertexEvent.cpp
 * @author Gernot Walzl
 * @date   2012-10-25
 */

#include "data/3d/skel/VertexEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelVertexData.h"

namespace data { namespace _3d { namespace skel {

VertexEvent::VertexEvent() {
    type_ = AbstractEvent::VERTEX_EVENT;
}

VertexEvent::~VertexEvent() {
    node_.reset();
    vertex_1_.reset();
    vertex_2_.reset();
    facet_1_.reset();
    facet_2_.reset();
}

VertexEventSPtr VertexEvent::create() {
    VertexEventSPtr result = VertexEventSPtr(new VertexEvent());
    return result;
}

NodeSPtr VertexEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void VertexEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double VertexEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

VertexSPtr VertexEvent::getVertex1() const {
    DEBUG_SPTR(vertex_1_);
    return vertex_1_;
}

void VertexEvent::setVertex1(VertexSPtr vertex_1) {
    this->vertex_1_ = vertex_1;
}

VertexSPtr VertexEvent::getVertex2() const {
    DEBUG_SPTR(vertex_2_);
    return vertex_2_;
}

void VertexEvent::setVertex2(VertexSPtr vertex_2) {
    this->vertex_2_ = vertex_2;
}

FacetSPtr VertexEvent::getFacet1() const {
    DEBUG_SPTR(facet_1_);
    return facet_1_;
}

void VertexEvent::setFacet1(FacetSPtr facet_1) {
    this->facet_1_ = facet_1;
}

FacetSPtr VertexEvent::getFacet2() const {
    DEBUG_SPTR(facet_2_);
    return facet_2_;
}

void VertexEvent::setFacet2(FacetSPtr facet_2) {
    this->facet_2_ = facet_2;
}

void VertexEvent::setHighlight(bool highlight) {
    if (!vertex_1_->hasData()) {
        SkelVertexData::create(vertex_1_);
    }
    vertex_1_->getData()->setHighlight(highlight);
    if (!vertex_2_->hasData()) {
        SkelVertexData::create(vertex_2_);
    }
    vertex_2_->getData()->setHighlight(highlight);
}

} } }
