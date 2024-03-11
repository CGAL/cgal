/**
 * @file   data/3d/skel/DblTriangleEvent.cpp
 * @author Gernot Walzl
 * @date   2012-09-11
 */

#include "data/3d/skel/DblTriangleEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

DblTriangleEvent::DblTriangleEvent() {
    type_ = AbstractEvent::DBL_TRIANGLE_EVENT;
}

DblTriangleEvent::~DblTriangleEvent() {
    node_.reset();
}

DblTriangleEventSPtr DblTriangleEvent::create() {
    DblTriangleEventSPtr result = DblTriangleEventSPtr(new DblTriangleEvent());
    return result;
}

NodeSPtr DblTriangleEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void DblTriangleEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double DblTriangleEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr DblTriangleEvent::getEdge() const {
    EdgeSPtr result = edge_;
    DEBUG_SPTR(result);
    return result;
}

void DblTriangleEvent::setEdge(EdgeSPtr edge) {
    this->edge_ = edge;
}

void DblTriangleEvent::getVertices(VertexSPtr out[4]) const {
    for (unsigned int i = 0; i < 4; i++) {
        out[i] = VertexSPtr();
    }
    FacetSPtr facet_l = edge_->getFacetL();
    FacetSPtr facet_r = edge_->getFacetR();
    out[0] = edge_->getVertexSrc();
    out[1] = edge_->getVertexDst();
    out[2] = edge_->next(facet_l)->dst(facet_l);
    out[3] = edge_->next(facet_r)->dst(facet_r);
}

void DblTriangleEvent::getEdges(EdgeSPtr out[5]) const {
    for (unsigned int i = 0; i < 5; i++) {
        out[i] = EdgeSPtr();
    }
    FacetSPtr facet_l = edge_->getFacetL();
    FacetSPtr facet_r = edge_->getFacetR();
    out[0] = edge_;
    out[1] = edge_->next(facet_l);
    out[2] = edge_->prev(facet_l);
    out[3] = edge_->next(facet_r);
    out[4] = edge_->prev(facet_r);
}

void DblTriangleEvent::setHighlight(bool highlight) {
    VertexSPtr vertices[4];
    getVertices(vertices);
    for (unsigned int i = 0; i < 4; i++) {
        if (!vertices[i]->hasData()) {
            SkelVertexData::create(vertices[i]);
        }
        vertices[i]->getData()->setHighlight(highlight);
    }
    EdgeSPtr edges[5];
    getEdges(edges);
    for (unsigned int i = 0; i < 5; i++) {
        if (!edges[i]->hasData()) {
            SkelEdgeData::create(edges[i]);
        }
        edges[i]->getData()->setHighlight(highlight);
    }
}

} } }
