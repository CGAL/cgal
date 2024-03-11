/**
 * @file   data/3d/skel/TetrahedronEvent.cpp
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#include "data/3d/skel/TetrahedronEvent.h"

#include "debug.h"
#include "data/3d/Edge.h"
#include "data/3d/Vertex.h"
#include "data/3d/Facet.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

TetrahedronEvent::TetrahedronEvent() {
    type_ = AbstractEvent::TETRAHEDRON_EVENT;
}

TetrahedronEvent::~TetrahedronEvent() {
    node_.reset();
}

TetrahedronEventSPtr TetrahedronEvent::create() {
    TetrahedronEventSPtr result = TetrahedronEventSPtr(new TetrahedronEvent());
    return result;
}

NodeSPtr TetrahedronEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void TetrahedronEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double TetrahedronEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr TetrahedronEvent::getEdgeBegin() const {
    EdgeSPtr result = edge_begin_;
    DEBUG_SPTR(result);
    return result;
}

void TetrahedronEvent::setEdgeBegin(EdgeSPtr edge_begin) {
    this->edge_begin_ = edge_begin;
}

void TetrahedronEvent::getVertices(VertexSPtr out[4]) const {
    for (unsigned int i = 0; i < 4; i++) {
        out[i] = VertexSPtr();
    }
    out[0] = edge_begin_->getVertexSrc();
    out[1] = edge_begin_->getVertexDst();
    EdgeSPtr edge_l = edge_begin_->next(edge_begin_->getFacetL());
    out[2] = edge_l->dst(edge_begin_->getFacetL());
    EdgeSPtr edge_r = edge_begin_->next(edge_begin_->getFacetR());
    out[3] = edge_r->dst(edge_begin_->getFacetR());
}

void TetrahedronEvent::getEdges(EdgeSPtr out[6]) const {
    for (unsigned int i = 0; i < 6; i++) {
        out[i] = EdgeSPtr();
    }
    out[0] = edge_begin_;
    out[1] = edge_begin_->prev(edge_begin_->getFacetL());
    out[2] = edge_begin_->next(edge_begin_->getFacetL());
    out[3] = edge_begin_->prev(edge_begin_->getFacetR());
    out[4] = edge_begin_->next(edge_begin_->getFacetR());
    FacetSPtr other = out[2]->other(edge_begin_->getFacetL());
    out[5] = out[2]->prev(other);
}

void TetrahedronEvent::getFacets(FacetSPtr out[4]) const {
    for (unsigned int i = 0; i < 4; i++) {
        out[i] = FacetSPtr();
    }
    out[0] = edge_begin_->getFacetL();
    out[1] = edge_begin_->getFacetR();
    out[2] = out[0]->prev(edge_begin_->getVertexDst());
    out[3] = out[1]->prev(edge_begin_->getVertexSrc());
}

void TetrahedronEvent::setHighlight(bool highlight) {
    VertexSPtr vertices[4];
    getVertices(vertices);
    for (unsigned int i = 0; i < 4; i++) {
        if (!vertices[i]->hasData()) {
            SkelVertexData::create(vertices[i]);
        }
        vertices[i]->getData()->setHighlight(highlight);
    }
//    EdgeSPtr edges[6];
//    getEdges(edges);
//    for (unsigned int i = 0; i < 6; i++) {
//        if (!edges[i]->hasData()) {
//            SkelEdgeData::create(edges[i]);
//        }
//        edges[i]->getData()->setHighlight(highlight);
//    }
}

} } }
