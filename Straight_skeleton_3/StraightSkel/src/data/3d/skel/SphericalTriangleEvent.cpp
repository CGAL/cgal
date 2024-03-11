/**
 * @file   data/3d/skel/SphericalTriangleEvent.cpp
 * @author Gernot Walzl
 * @date   2012-11-29
 */

#include "data/3d/skel/SphericalTriangleEvent.h"

#include "debug.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/SphericalSkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

SphericalTriangleEvent::SphericalTriangleEvent() {
    type_ = SphericalAbstractEvent::TRIANGLE_EVENT;
}

SphericalTriangleEvent::~SphericalTriangleEvent() {
    node_.reset();
    edge_begin_.reset();
}

SphericalTriangleEventSPtr SphericalTriangleEvent::create() {
    SphericalTriangleEventSPtr result = SphericalTriangleEventSPtr(new SphericalTriangleEvent());
    return result;
}

CircularNodeSPtr SphericalTriangleEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void SphericalTriangleEvent::setNode(CircularNodeSPtr node) {
    this->node_ = node;
}

double SphericalTriangleEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

CircularEdgeSPtr SphericalTriangleEvent::getEdgeBegin() const {
    DEBUG_SPTR(edge_begin_);
    return edge_begin_;
}

void SphericalTriangleEvent::setEdgeBegin(CircularEdgeSPtr edge_begin) {
    this->edge_begin_ = edge_begin;
}

void SphericalTriangleEvent::getVertices(CircularVertexSPtr out[3]) const {
    for (unsigned int i = 0; i < 3; i++) {
        out[i] = CircularVertexSPtr();
    }
    out[0] = edge_begin_->getVertexSrc();
    out[1] = edge_begin_->getVertexDst();
    out[2] = edge_begin_->next()->getVertexDst();
}

void SphericalTriangleEvent::getEdges(CircularEdgeSPtr out[3]) const {
    for (unsigned int i = 0; i < 3; i++) {
        out[i] = CircularEdgeSPtr();
    }
    out[0] = edge_begin_;
    out[1] = edge_begin_->next();
    out[2] = edge_begin_->prev();
}

void SphericalTriangleEvent::setHighlight(bool highlight) {
    CircularEdgeSPtr edges[3];
    getEdges(edges);
    for (unsigned int i = 0; i < 3; i++) {
        if (!edges[i]->hasData()) {
            SphericalSkelEdgeData::create(edges[i]);
        }
        edges[i]->getData()->setHighlight(highlight);
    }
}

} } }
