/**
 * @file   data/3d/skel/DblEdgeMergeEvent.cpp
 * @author Gernot Walzl
 * @date   2012-10-30
 */

#include "data/3d/skel/DblEdgeMergeEvent.h"

#include "debug.h"
#include "data/3d/Edge.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

DblEdgeMergeEvent::DblEdgeMergeEvent() {
    type_ = AbstractEvent::DBL_EDGE_MERGE_EVENT;
}

DblEdgeMergeEvent::~DblEdgeMergeEvent() {
    node_.reset();
    facet_1_.reset();
    edge_11_.reset();
    edge_12_.reset();
    facet_2_.reset();
    edge_21_.reset();
    edge_22_.reset();
}

DblEdgeMergeEventSPtr DblEdgeMergeEvent::create() {
    DblEdgeMergeEventSPtr result = DblEdgeMergeEventSPtr(new DblEdgeMergeEvent());
    return result;
}

NodeSPtr DblEdgeMergeEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void DblEdgeMergeEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double DblEdgeMergeEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

FacetSPtr DblEdgeMergeEvent::getFacet1() const {
    DEBUG_SPTR(facet_1_);
    return facet_1_;
}

void DblEdgeMergeEvent::setFacet1(FacetSPtr facet_1) {
    this->facet_1_ = facet_1;
}

EdgeSPtr DblEdgeMergeEvent::getEdge11() const {
    DEBUG_SPTR(edge_11_);
    return edge_11_;
}

void DblEdgeMergeEvent::setEdge11(EdgeSPtr edge_11) {
    this->edge_11_ = edge_11;
}

EdgeSPtr DblEdgeMergeEvent::getEdge12() const {
    DEBUG_SPTR(edge_12_);
    return edge_12_;
}

void DblEdgeMergeEvent::setEdge12(EdgeSPtr edge_12) {
    this->edge_12_ = edge_12;
}

FacetSPtr DblEdgeMergeEvent::getFacet2() const {
    DEBUG_SPTR(facet_2_);
    return facet_2_;
}

void DblEdgeMergeEvent::setFacet2(FacetSPtr facet_2) {
    this->facet_2_ = facet_2;
}

EdgeSPtr DblEdgeMergeEvent::getEdge21() const {
    DEBUG_SPTR(edge_21_);
    return edge_21_;
}

void DblEdgeMergeEvent::setEdge21(EdgeSPtr edge_21) {
    this->edge_21_ = edge_21;
}

EdgeSPtr DblEdgeMergeEvent::getEdge22() const {
    DEBUG_SPTR(edge_22_);
    return edge_22_;
}

void DblEdgeMergeEvent::setEdge22(EdgeSPtr edge_22) {
    this->edge_22_ = edge_22;
}

void DblEdgeMergeEvent::getVertices(VertexSPtr out[4]) const {
    for (unsigned int i = 0; i < 4; i++) {
        out[i] = VertexSPtr();
    }
    out[0] = edge_11_->dst(facet_1_);
    out[1] = edge_21_->dst(facet_2_);
    out[2] = edge_12_->src(facet_1_);
    out[3] = edge_22_->src(facet_2_);
}

void DblEdgeMergeEvent::getEdges(EdgeSPtr out[4]) const {
    for (unsigned int i = 0; i < 4; i++) {
        out[i] = EdgeSPtr();
    }
    out[0] = edge_11_->next(facet_1_);
    out[1] = out[0]->next(facet_1_);
    FacetSPtr facet_other = edge_11_->other(facet_1_);
    out[2] = edge_12_->next(facet_other);
    out[3] = out[2]->next(facet_other);
}

void DblEdgeMergeEvent::setHighlight(bool highlight) {
    if (!edge_11_->hasData()) {
        SkelEdgeData::create(edge_11_);
    }
    edge_11_->getData()->setHighlight(highlight);
    if (!edge_12_->hasData()) {
        SkelEdgeData::create(edge_12_);
    }
    edge_12_->getData()->setHighlight(highlight);
    if (!edge_21_->hasData()) {
        SkelEdgeData::create(edge_21_);
    }
    edge_21_->getData()->setHighlight(highlight);
    if (!edge_22_->hasData()) {
        SkelEdgeData::create(edge_22_);
    }
    edge_22_->getData()->setHighlight(highlight);
    EdgeSPtr edges[4];
    getEdges(edges);
    for (unsigned int i = 0; i < 4; i++) {
        if (!edges[i]->hasData()) {
            SkelEdgeData::create(edges[i]);
        }
        edges[i]->getData()->setHighlight(highlight);
    }
}

} } }
