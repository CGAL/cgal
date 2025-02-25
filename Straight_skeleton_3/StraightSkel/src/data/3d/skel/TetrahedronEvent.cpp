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
#include "util/StringFactory.h"

#include <sstream>
#include <string>

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

CGAL::FT TetrahedronEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr TetrahedronEvent::getEdgeBegin() const {
    DEBUG_WPTR(edge_begin_);
    return edge_begin_.lock();
}

void TetrahedronEvent::setEdgeBegin(EdgeSPtr edge_begin) {
    this->edge_begin_ = edge_begin;
    this->neighborhood_ = EdgeFacetNeighborhood(edge_begin);
}

void TetrahedronEvent::getVertices(VertexSPtr out[4]) const {
    EdgeSPtr edge_begin = getEdgeBegin();

    for (unsigned int i = 0; i < 4; i++) {
        out[i] = VertexSPtr();
    }
    out[0] = edge_begin->getVertexSrc();
    out[1] = edge_begin->getVertexDst();
    EdgeSPtr edge_l = edge_begin->next(edge_begin->getFacetL());
    out[2] = edge_l->dst(edge_begin->getFacetL());
    EdgeSPtr edge_r = edge_begin->next(edge_begin->getFacetR());
    out[3] = edge_r->dst(edge_begin->getFacetR());
}

void TetrahedronEvent::getEdges(EdgeSPtr out[6]) const {
    EdgeSPtr edge_begin = getEdgeBegin();

    for (unsigned int i = 0; i < 6; i++) {
        out[i] = EdgeSPtr();
    }
    out[0] = edge_begin;
    out[1] = edge_begin->prev(edge_begin->getFacetL());
    out[2] = edge_begin->next(edge_begin->getFacetL());
    out[3] = edge_begin->prev(edge_begin->getFacetR());
    out[4] = edge_begin->next(edge_begin->getFacetR());
    FacetSPtr other = out[2]->other(edge_begin->getFacetL());
    out[5] = out[2]->prev(other);
}

void TetrahedronEvent::getFacets(FacetSPtr out[4]) const {
    EdgeSPtr edge_begin = getEdgeBegin();

    for (unsigned int i = 0; i < 4; i++) {
        out[i] = FacetSPtr();
    }
    out[0] = edge_begin->getFacetL();
    out[1] = edge_begin->getFacetR();
    out[2] = out[0]->prev(edge_begin->getVertexDst());
    out[3] = out[1]->prev(edge_begin->getVertexSrc());
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

std::string TetrahedronEvent::toString() const {
    VertexSPtr vertices[4];
    getVertices(vertices);

    EdgeSPtr edges[6];
    getEdges(edges);

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "TetrahedronEvent\n";
    sstr << "\t(ID=" << getID() << "; step ID=" << getStepID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(vertices";
    for(int i=0; i<4; ++i)
        sstr << " " << vertices[i]->getID();
    sstr << ")\n";
    sstr << "\t(edges";
    for(int i=0; i<6; ++i)
        sstr << " " << edges[i]->getID();
    sstr << ")";
    return sstr.str();
}

bool TetrahedronEvent::isValid() const {
    return node_ && !edge_begin_.expired();
}

bool TetrahedronEvent::isObsolete() const {
    if (EdgeSPtr edge = getEdgeBegin()) {
        return ! neighborhood_.checkNeighborhoodConsistency(edge);
    }

    return false;
}

} } }
