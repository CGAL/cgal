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
#include "util/StringFactory.h"

#include <sstream>
#include <string>

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
    CGAL_SS3_DEBUG_SPTR(node_);
    return node_;
}

void TriangleEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT TriangleEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

FacetSPtr TriangleEvent::getFacet() const {
    CGAL_SS3_DEBUG_WPTR(facet_);
    return facet_.lock();
}

void TriangleEvent::setFacet(FacetSPtr facet) {
    this->facet_ = facet;
}

EdgeSPtr TriangleEvent::getEdgeBegin() const {
    CGAL_SS3_DEBUG_WPTR(edge_begin_);
    return edge_begin_.lock();
}

void TriangleEvent::setEdgeBegin(EdgeSPtr edge_begin) {
    this->edge_begin_ = edge_begin;
    this->neighborhood_ = EdgeFacetNeighborhood(edge_begin);
}

void TriangleEvent::getVertices(VertexSPtr out[3]) const {
    FacetSPtr facet = getFacet();
    EdgeSPtr edge_begin = getEdgeBegin();

    for (unsigned int i = 0; i < 3; i++) {
        out[i] = VertexSPtr();
    }
    out[0] = edge_begin->src(facet);
    out[1] = edge_begin->dst(facet);
    out[2] = edge_begin->next(facet)->dst(facet);
}

void TriangleEvent::getEdges(EdgeSPtr out[3]) const {
    FacetSPtr facet = getFacet();
    EdgeSPtr edge_begin = getEdgeBegin();

    for (unsigned int i = 0; i < 3; i++) {
        out[i] = EdgeSPtr();
    }
    out[0] = edge_begin;
    out[1] = edge_begin->next(facet);
    out[2] = edge_begin->prev(facet);
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

std::string TriangleEvent::toString() const {
    FacetSPtr facet = getFacet();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "TriangleEvent\n";
    sstr << "\t(ID=" << getID() << "; step ID=" << getStepID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(facet=" << facet->getID() << ")";
    return sstr.str();
}

bool TriangleEvent::isValid() const {
    return node_ && !facet_.expired() && !edge_begin_.expired();
}

bool TriangleEvent::isObsolete() const {
    if (EdgeSPtr edge = getEdgeBegin()) {
        return ! neighborhood_.checkNeighborhoodConsistency(edge);
    }

    return false;
}

bool TriangleEvent::operator==(const TriangleEvent& other) const {
    return (node_->getOffset() == other.node_->getOffset()) &&
           (*(node_->getPoint()) == *(other.node_->getPoint())) &&
           (facet_.lock() == other.facet_.lock()) &&
           (edge_begin_.lock() == other.edge_begin_.lock());
}

} } }
