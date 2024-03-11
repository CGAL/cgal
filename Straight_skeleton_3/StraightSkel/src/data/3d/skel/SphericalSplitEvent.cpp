/**
 * @file   data/3d/skel/SphericalSplitEvent.cpp
 * @author Gernot Walzl
 * @date   2012-11-29
 */

#include "data/3d/skel/SphericalSplitEvent.h"

#include "debug.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/SphericalSkelVertexData.h"
#include "data/3d/skel/SphericalSkelEdgeData.h"

namespace data { namespace _3d { namespace skel {

SphericalSplitEvent::SphericalSplitEvent() {
    type_ = SphericalAbstractEvent::SPLIT_EVENT;
}

SphericalSplitEvent::~SphericalSplitEvent() {
    node_.reset();
    vertex_.reset();
    edge_.reset();
}

SphericalSplitEventSPtr SphericalSplitEvent::create() {
    SphericalSplitEventSPtr result = SphericalSplitEventSPtr(new SphericalSplitEvent());
    return result;
}

CircularNodeSPtr SphericalSplitEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void SphericalSplitEvent::setNode(CircularNodeSPtr node) {
    this->node_ = node;
}

double SphericalSplitEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

CircularVertexSPtr SphericalSplitEvent::getVertex() const {
    DEBUG_SPTR(vertex_);
    return vertex_;
}

void SphericalSplitEvent::setVertex(CircularVertexSPtr vertex) {
    this->vertex_ = vertex;
}

CircularEdgeSPtr SphericalSplitEvent::getEdge() const {
    DEBUG_SPTR(edge_);
    return edge_;
}

void SphericalSplitEvent::setEdge(CircularEdgeSPtr edge) {
    this->edge_ = edge;
}

void SphericalSplitEvent::setHighlight(bool highlight) {
    if (!vertex_->hasData()) {
        SphericalSkelVertexData::create(vertex_);
    }
    vertex_->getData()->setHighlight(highlight);
    if (!edge_->hasData()) {
        SphericalSkelEdgeData::create(edge_);
    }
    edge_->getData()->setHighlight(highlight);
}

} } }
