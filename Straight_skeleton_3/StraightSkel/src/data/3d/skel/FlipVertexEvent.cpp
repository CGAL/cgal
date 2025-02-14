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
 * @file   data/3d/skel/FlipVertexEvent.cpp
 * @author Gernot Walzl
 * @date   2012-11-22
 */

#include "data/3d/skel/FlipVertexEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Facet.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelVertexData.h"
#include "util/StringFactory.h"

#include <sstream>

namespace data { namespace _3d { namespace skel {

FlipVertexEvent::FlipVertexEvent() {
    type_ = AbstractEvent::FLIP_VERTEX_EVENT;
}

FlipVertexEvent::~FlipVertexEvent() {
    node_.reset();
    vertex_1_.reset();
    vertex_2_.reset();
    facet_1_.reset();
    facet_2_.reset();
}

FlipVertexEventSPtr FlipVertexEvent::create() {
    FlipVertexEventSPtr result = FlipVertexEventSPtr(new FlipVertexEvent());
    return result;
}

NodeSPtr FlipVertexEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void FlipVertexEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT FlipVertexEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

VertexSPtr FlipVertexEvent::getVertex1() const {
    DEBUG_WPTR(vertex_1_);
    return vertex_1_.lock();
}

void FlipVertexEvent::setVertex1(VertexSPtr vertex_1) {
    this->vertex_1_ = vertex_1;
}

VertexSPtr FlipVertexEvent::getVertex2() const {
    DEBUG_WPTR(vertex_2_);
    return vertex_2_.lock();
}

void FlipVertexEvent::setVertex2(VertexSPtr vertex_2) {
    this->vertex_2_ = vertex_2;
}

FacetSPtr FlipVertexEvent::getFacet1() const {
    DEBUG_WPTR(facet_1_);
    return facet_1_.lock();
}

void FlipVertexEvent::setFacet1(FacetSPtr facet_1) {
    this->facet_1_ = facet_1;
}

FacetSPtr FlipVertexEvent::getFacet2() const {
    DEBUG_WPTR(facet_2_);
    return facet_2_.lock();
}

void FlipVertexEvent::setFacet2(FacetSPtr facet_2) {
    this->facet_2_ = facet_2;
}

void FlipVertexEvent::setHighlight(bool highlight) {
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

std::string FlipVertexEvent::toString() const {
    VertexSPtr vertex_1 = getVertex1();
    VertexSPtr vertex_2 = getVertex2();
    FacetSPtr facet_1 = getFacet1();
    FacetSPtr facet_2 = getFacet2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "FlipVertexEvent\n";
    sstr << "\t(ID=" << getID() << "; step ID=" << getStepID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(vertex A=" << vertex_1->getID() << "; vertex B=" << vertex_2->getID() << ")\n";
    sstr << "\t(facet A=" << facet_1->getID() << "; facet B=" << facet_2->getID() << ")";
    return sstr.str();
}

bool FlipVertexEvent::isValid() const {
    return node_ && !vertex_1_.expired() && !vertex_2_.expired() && !facet_1_.expired() && !facet_2_.expired();
}

} } }
