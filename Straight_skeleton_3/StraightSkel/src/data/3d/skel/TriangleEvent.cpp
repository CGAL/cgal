/**
 * @file   data/3d/skel/TriangleEvent.cpp
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#include "data/3d/skel/TriangleEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

TriangleEvent::TriangleEvent() {
    type_ = AbstractEvent::TRIANGLE_EVENT;
}

TriangleEvent::~TriangleEvent() {
    node_.reset();
}

TriangleEventSPtr TriangleEvent::create() {
    TriangleEventSPtr result = TriangleEventSPtr(new TriangleEvent());
    return result;
}

NodeSPtr TriangleEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void TriangleEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double TriangleEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

FacetSPtr TriangleEvent::getFacet() const {
    FacetSPtr result = facet_;
    DEBUG_SPTR(result);
    return result;
}

void TriangleEvent::setFacet(FacetSPtr facet) {
    this->facet_ = facet;
}

EdgeSPtr TriangleEvent::getEdgeBegin() const {
    EdgeSPtr result = edge_begin_;
    DEBUG_SPTR(result);
    return result;
}

void TriangleEvent::setEdgeBegin(EdgeSPtr edge_begin) {
    this->edge_begin_ = edge_begin;
}

void TriangleEvent::getVertices(VertexSPtr out[3]) const {
    for (unsigned int i = 0; i < 3; i++) {
        out[i] = VertexSPtr();
    }
    out[0] = edge_begin_->src(facet_);
    out[1] = edge_begin_->dst(facet_);
    out[2] = edge_begin_->next(facet_)->dst(facet_);
}

void TriangleEvent::getEdges(EdgeSPtr out[3]) const {
    for (unsigned int i = 0; i < 3; i++) {
        out[i] = EdgeSPtr();
    }
    out[0] = edge_begin_;
    out[1] = edge_begin_->next(facet_);
    out[2] = edge_begin_->prev(facet_);
}

void TriangleEvent::setHighlight(bool highlight) {
    VertexSPtr vertices[3];
    getVertices(vertices);
    for (unsigned int i = 0; i < 3; i++) {
        if (!vertices[i]->hasData()) {
            SkelVertexData::create(vertices[i]);
        }
        vertices[i]->getData()->setHighlight(highlight);
    }
    EdgeSPtr edges[3];
    getEdges(edges);
    for (unsigned int i = 0; i < 3; i++) {
        if (!edges[i]->hasData()) {
            SkelEdgeData::create(edges[i]);
        }
        edges[i]->getData()->setHighlight(highlight);
    }
}

} } }
