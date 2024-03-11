/**
 * @file   data/3d/skel/SplitMergeEvent.cpp
 * @author Gernot Walzl
 * @date   2013-08-09
 */

#include "data/3d/skel/SplitMergeEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelVertexData.h"

namespace data { namespace _3d { namespace skel {

SplitMergeEvent::SplitMergeEvent() {
    type_ = AbstractEvent::SPLIT_MERGE_EVENT;
}

SplitMergeEvent::~SplitMergeEvent() {
    node_.reset();
    vertex_1_.reset();
    vertex_2_.reset();
    facet_1_.reset();
    facet_2_.reset();
}

SplitMergeEventSPtr SplitMergeEvent::create() {
    SplitMergeEventSPtr result = SplitMergeEventSPtr(new SplitMergeEvent());
    return result;
}

NodeSPtr SplitMergeEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void SplitMergeEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double SplitMergeEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

VertexSPtr SplitMergeEvent::getVertex1() const {
    DEBUG_SPTR(vertex_1_);
    return vertex_1_;
}

void SplitMergeEvent::setVertex1(VertexSPtr vertex_1) {
    this->vertex_1_ = vertex_1;
}

VertexSPtr SplitMergeEvent::getVertex2() const {
    DEBUG_SPTR(vertex_2_);
    return vertex_2_;
}

void SplitMergeEvent::setVertex2(VertexSPtr vertex_2) {
    this->vertex_2_ = vertex_2;
}

FacetSPtr SplitMergeEvent::getFacet1() const {
    DEBUG_SPTR(facet_1_);
    return facet_1_;
}

void SplitMergeEvent::setFacet1(FacetSPtr facet_1) {
    this->facet_1_ = facet_1;
}

FacetSPtr SplitMergeEvent::getFacet2() const {
    DEBUG_SPTR(facet_2_);
    return facet_2_;
}

void SplitMergeEvent::setFacet2(FacetSPtr facet_2) {
    this->facet_2_ = facet_2;
}

void SplitMergeEvent::setHighlight(bool highlight) {
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
