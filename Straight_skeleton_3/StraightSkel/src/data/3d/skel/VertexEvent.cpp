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
 * @file   data/3d/skel/VertexEvent.cpp
 * @author Gernot Walzl
 * @date   2012-10-25
 */

#include "data/3d/skel/VertexEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Facet.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelVertexData.h"
#include "util/StringFactory.h"

#include <sstream>

namespace data { namespace _3d { namespace skel {

VertexEvent::VertexEvent() {
    type_ = AbstractEvent::VERTEX_EVENT;
}

VertexEvent::~VertexEvent() {
    node_.reset();
    vertex_1_.reset();
    vertex_2_.reset();
    facet_1_.reset();
    facet_2_.reset();
}

VertexEventSPtr VertexEvent::create() {
    VertexEventSPtr result = VertexEventSPtr(new VertexEvent());
    return result;
}

NodeSPtr VertexEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void VertexEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT VertexEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

VertexSPtr VertexEvent::getVertex1() const {
    DEBUG_WPTR(vertex_1_);
    return vertex_1_.lock();
}

void VertexEvent::setVertex1(VertexSPtr vertex_1) {
    this->vertex_1_ = vertex_1;
}

VertexSPtr VertexEvent::getVertex2() const {
    DEBUG_WPTR(vertex_2_);
    return vertex_2_.lock();
}

void VertexEvent::setVertex2(VertexSPtr vertex_2) {
    this->vertex_2_ = vertex_2;
}

FacetSPtr VertexEvent::getFacet1() const {
    DEBUG_WPTR(facet_1_);
    return facet_1_.lock();
}

void VertexEvent::setFacet1(FacetSPtr facet_1) {
    this->facet_1_ = facet_1;
}

FacetSPtr VertexEvent::getFacet2() const {
    DEBUG_WPTR(facet_2_);
    return facet_2_.lock();
}

void VertexEvent::setFacet2(FacetSPtr facet_2) {
    this->facet_2_ = facet_2;
}

void VertexEvent::setHighlight(bool highlight) {
    VertexSPtr vertex_1 = getVertex1();
    VertexSPtr vertex_2 = getVertex2();

    if (!vertex_1->hasData()) {
        SkelVertexData::create(vertex_1);
    }
    vertex_1->getData()->setHighlight(highlight);
    if (!vertex_2->hasData()) {
        SkelVertexData::create(vertex_2);
    }
    vertex_2->getData()->setHighlight(highlight);
}

std::string VertexEvent::toString() const {
    VertexSPtr vertex_1 = getVertex1();
    VertexSPtr vertex_2 = getVertex2();
    FacetSPtr facet_1 = getFacet1();
    FacetSPtr facet_2 = getFacet2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "VertexEvent\n";
    sstr << "\t(ID=" << getID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(facet A=" << facet_1->getID()
         << "; vertex A=" << vertex_1->getID() << ")\n";
    sstr << "\t(facet B=" << facet_2->getID()
         << "; vertex B2=" << vertex_2->getID() << ")";
    return sstr.str();
}

bool VertexEvent::isValid() const {
    return node_ && !vertex_1_.expired() && !vertex_2_.expired() && !facet_1_.expired() && !facet_2_.expired();
}

} } }
