/**
 * @file   data/3d/skel/PierceEvent.cpp
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#include "data/3d/skel/PierceEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "data/3d/skel/SkelFacetData.h"
#include <list>

namespace data { namespace _3d { namespace skel {

PierceEvent::PierceEvent() {
    type_ = AbstractEvent::PIERCE_EVENT;
}

PierceEvent::~PierceEvent() {
    node_.reset();
    facet_.reset();
}

PierceEventSPtr PierceEvent::create() {
    PierceEventSPtr result = PierceEventSPtr(new PierceEvent());
    return result;
}

NodeSPtr PierceEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void PierceEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

double PierceEvent::getOffset() const {
    double result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

FacetSPtr PierceEvent::getFacet() const {
    DEBUG_SPTR(facet_);
    return facet_;
}

void PierceEvent::setFacet(FacetSPtr facet) {
    this->facet_ = facet;
}

VertexSPtr PierceEvent::getVertex() const {
    DEBUG_SPTR(vertex_);
    return vertex_;
}

void PierceEvent::setVertex(VertexSPtr vertex) {
    this->vertex_ = vertex;
}

void PierceEvent::setHighlight(bool highlight) {
    if (!vertex_->hasData()) {
        SkelVertexData::create(vertex_);
    }
    vertex_->getData()->setHighlight(highlight);
    if (!facet_->hasData()) {
        SkelFacetData::create(facet_);
    }
    facet_->getData()->setHighlight(highlight);
    std::list<EdgeSPtr>::iterator it_e = facet_->edges().begin();
    while (it_e != facet_->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (!edge->hasData()) {
            SkelEdgeData::create(edge);
        }
        edge->getData()->setHighlight(highlight);
    }
}

} } }
