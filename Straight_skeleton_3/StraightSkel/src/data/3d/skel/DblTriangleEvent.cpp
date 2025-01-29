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
#include "util/StringFactory.h"

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

CGAL::FT DblTriangleEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr DblTriangleEvent::getEdge() const {
    DEBUG_WPTR(edge_);
    return edge_.lock();
}

void DblTriangleEvent::setEdge(EdgeSPtr edge) {
    this->edge_ = edge;
}

void DblTriangleEvent::getVertices(VertexSPtr out[4]) const {
    EdgeSPtr edge = getEdge();

    for (unsigned int i = 0; i < 4; i++) {
        out[i] = VertexSPtr();
    }
    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    out[0] = edge->getVertexSrc();
    out[1] = edge->getVertexDst();
    out[2] = edge->next(facet_l)->dst(facet_l);
    out[3] = edge->next(facet_r)->dst(facet_r);
}

void DblTriangleEvent::getEdges(EdgeSPtr out[5]) const {
    EdgeSPtr edge = getEdge();

    for (unsigned int i = 0; i < 5; i++) {
        out[i] = EdgeSPtr();
    }
    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    out[0] = edge;
    out[1] = edge->next(facet_l);
    out[2] = edge->prev(facet_l);
    out[3] = edge->next(facet_r);
    out[4] = edge->prev(facet_r);
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

std::string DblTriangleEvent::toString() const {
    EdgeSPtr edge = getEdge();

    VertexSPtr vertices[4];
    getVertices(vertices);

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "DblTriangleEvent\n";
    sstr << "\t(ID=" << getID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(edge=" << edge->getID() << "\n\t\t[" << edge->getVertexSrc()->toString() << "\n\t\t "
                                                     << edge->getVertexDst()->toString() << "])";
    return sstr.str();
}

bool DblTriangleEvent::isValid() const {
    return node_ && !edge_.expired();
}

} } }
