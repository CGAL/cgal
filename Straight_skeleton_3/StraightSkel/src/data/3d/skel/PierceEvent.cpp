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
#include "util/StringFactory.h"

#include <list>
#include <sstream>

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
    CGAL_SS3_DEBUG_SPTR(node_);
    return node_;
}

void PierceEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT PierceEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

FacetSPtr PierceEvent::getFacet() const {
    CGAL_SS3_DEBUG_WPTR(facet_);
    return facet_.lock();
}

void PierceEvent::setFacet(FacetSPtr facet) {
    this->facet_ = facet;
}

VertexSPtr PierceEvent::getVertex() const {
    CGAL_SS3_DEBUG_WPTR(vertex_);
    return vertex_.lock();
}

void PierceEvent::setVertex(VertexSPtr vertex) {
    this->vertex_ = vertex;
    this->neighborhood_ = VertexFacetNeighborhood(vertex);
}

void PierceEvent::setHighlight(bool highlight) {
    FacetSPtr facet = getFacet();
    VertexSPtr vertex = getVertex();

    if (!vertex->hasData()) {
        SkelVertexData::create(vertex);
    }
    vertex->getData()->setHighlight(highlight);
    if (!facet->hasData()) {
        SkelFacetData::create(facet);
    }
    facet->getData()->setHighlight(highlight);
    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (!edge->hasData()) {
            SkelEdgeData::create(edge);
        }
        edge->getData()->setHighlight(highlight);
    }
}

std::string PierceEvent::toString() const {
    FacetSPtr facet = getFacet();
    VertexSPtr vertex = getVertex();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "PierceEvent\n";
    sstr << "\t(ID=" << getID() << "; step ID=" << getStepID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(vertex=" << vertex->toString() << ")\n";
    sstr << "\t(face=" << facet->getID() << ")";
    return sstr.str();
}

bool PierceEvent::isValid() const {
    return node_ && !facet_.expired() && !vertex_.expired();
}

bool PierceEvent::isObsolete() const {
  if (VertexSPtr vertex = getVertex()) {
      // std::cout << "isObsolete(v" << vertex->getID() << ")?" << std::endl;
      return ! neighborhood_.checkNeighborhoodConsistency(vertex);
  }

  return false;
}

bool PierceEvent::operator==(const PierceEvent& other) const {
    return (node_->getOffset() == other.node_->getOffset()) &&
           (*(node_->getPoint()) == *(other.node_->getPoint())) &&
           (facet_.lock() == other.facet_.lock()) &&
           (vertex_.lock() == other.vertex_.lock());
}

} } }
