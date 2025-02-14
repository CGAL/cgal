// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   data/3d/skel/DblEdgeMergeEvent.cpp
 * @author Gernot Walzl
 * @date   2012-10-30
 */

#include "data/3d/skel/DblEdgeMergeEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "util/StringFactory.h"

#include <sstream>

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

CGAL::FT DblEdgeMergeEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

FacetSPtr DblEdgeMergeEvent::getFacet1() const {
    DEBUG_WPTR(facet_1_);
    return facet_1_.lock();
}

void DblEdgeMergeEvent::setFacet1(FacetSPtr facet_1) {
    this->facet_1_ = facet_1;
}

EdgeSPtr DblEdgeMergeEvent::getEdge11() const {
    DEBUG_WPTR(edge_11_);
    return edge_11_.lock();
}

void DblEdgeMergeEvent::setEdge11(EdgeSPtr edge_11) {
    this->edge_11_ = edge_11;
}

EdgeSPtr DblEdgeMergeEvent::getEdge12() const {
    DEBUG_WPTR(edge_12_);
    return edge_12_.lock();
}

void DblEdgeMergeEvent::setEdge12(EdgeSPtr edge_12) {
    this->edge_12_ = edge_12;
}

FacetSPtr DblEdgeMergeEvent::getFacet2() const {
    DEBUG_WPTR(facet_2_);
    return facet_2_.lock();
}

void DblEdgeMergeEvent::setFacet2(FacetSPtr facet_2) {
    this->facet_2_ = facet_2;
}

EdgeSPtr DblEdgeMergeEvent::getEdge21() const {
    DEBUG_WPTR(edge_21_);
    return edge_21_.lock();
}

void DblEdgeMergeEvent::setEdge21(EdgeSPtr edge_21) {
    this->edge_21_ = edge_21;
}

EdgeSPtr DblEdgeMergeEvent::getEdge22() const {
    DEBUG_WPTR(edge_22_);
    return edge_22_.lock();
}

void DblEdgeMergeEvent::setEdge22(EdgeSPtr edge_22) {
    this->edge_22_ = edge_22;
}

void DblEdgeMergeEvent::getVertices(VertexSPtr out[4]) const {
    EdgeSPtr edge_11 = getEdge11();
    EdgeSPtr edge_21 = getEdge21();
    EdgeSPtr edge_12 = getEdge12();
    EdgeSPtr edge_22 = getEdge22();
    FacetSPtr facet_1 = getFacet1();
    FacetSPtr facet_2 = getFacet2();

    for (unsigned int i = 0; i < 4; i++) {
        out[i] = VertexSPtr();
    }
    out[0] = edge_11->dst(facet_1);
    out[1] = edge_21->dst(facet_2);
    out[2] = edge_12->src(facet_1);
    out[3] = edge_22->src(facet_2);
}

void DblEdgeMergeEvent::getEdges(EdgeSPtr out[4]) const {
    EdgeSPtr edge_11 = getEdge11();
    EdgeSPtr edge_12 = getEdge12();
    FacetSPtr facet_1 = getFacet1();

    for (unsigned int i = 0; i < 4; i++) {
        out[i] = EdgeSPtr();
    }
    out[0] = edge_11->next(facet_1);
    out[1] = out[0]->next(facet_1);
    FacetSPtr facet_other = edge_11->other(facet_1);
    out[2] = edge_12->next(facet_other);
    out[3] = out[2]->next(facet_other);
}

void DblEdgeMergeEvent::setHighlight(bool highlight) {
    EdgeSPtr edge_11 = getEdge11();
    EdgeSPtr edge_21 = getEdge21();
    EdgeSPtr edge_12 = getEdge12();
    EdgeSPtr edge_22 = getEdge22();

    if (!edge_11->hasData()) {
        SkelEdgeData::create(edge_11);
    }
    edge_11->getData()->setHighlight(highlight);
    if (!edge_12->hasData()) {
        SkelEdgeData::create(edge_12);
    }
    edge_12->getData()->setHighlight(highlight);
    if (!edge_21->hasData()) {
        SkelEdgeData::create(edge_21);
    }
    edge_21->getData()->setHighlight(highlight);
    if (!edge_22->hasData()) {
        SkelEdgeData::create(edge_22);
    }
    edge_22->getData()->setHighlight(highlight);
    EdgeSPtr edges[4];
    getEdges(edges);
    for (unsigned int i = 0; i < 4; i++) {
        if (!edges[i]->hasData()) {
            SkelEdgeData::create(edges[i]);
        }
        edges[i]->getData()->setHighlight(highlight);
    }
}

std::string DblEdgeMergeEvent::toString() const {
    VertexSPtr vertices[4];
    getVertices(vertices);

    EdgeSPtr edges[4];
    getEdges(edges);

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "DblEdgeMergeEvent\n";
    sstr << "\t(ID=" << getID() << "; step ID=" << getStepID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    for (unsigned int i = 0; i < 4; i++) {
        sstr << "\t(vertex" << i << "=" << *(vertices[i]->getPoint()) << ")\n";
    }
    for (unsigned int i = 0; i < 4; i++) {
        sstr << "\t(edge" << i << "=" << edges[i]->getID() << ")\n";
    }

    return sstr.str();
}

bool DblEdgeMergeEvent::isValid() const {
    return node_ && !facet_1_.expired() && !edge_11_.expired() && !edge_12_.expired() && !facet_2_.expired() && !edge_21_.expired() && !edge_22_.expired();
}

} } }
