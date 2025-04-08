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
 * @file   data/3d/skel/SplitMergeEvent.cpp
 * @author Gernot Walzl
 * @date   2013-08-09
 */

#include "data/3d/skel/SplitMergeEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Facet.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelVertexData.h"
#include "util/StringFactory.h"

namespace data { namespace _3d { namespace skel {

SplitMergeEvent::SplitMergeEvent() {
    type_ = AbstractEvent::SPLIT_MERGE_EVENT;
}

SplitMergeEvent::~SplitMergeEvent() {
    node_.reset();
    vertex_1_.reset();
    vertex_2_.reset();
    facet_1_.reset();
    facet_2_.reset();
}

SplitMergeEventSPtr SplitMergeEvent::create() {
    SplitMergeEventSPtr result = SplitMergeEventSPtr(new SplitMergeEvent());
    return result;
}

NodeSPtr SplitMergeEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void SplitMergeEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT SplitMergeEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

VertexSPtr SplitMergeEvent::getVertex1() const {
    DEBUG_WPTR(vertex_1_);
    return vertex_1_.lock();
}

void SplitMergeEvent::setVertex1(VertexSPtr vertex_1) {
    this->vertex_1_ = vertex_1;
    this->neighborhood_1_ = VertexFacetNeighborhood(vertex_1);
}

VertexSPtr SplitMergeEvent::getVertex2() const {
    DEBUG_WPTR(vertex_2_);
    return vertex_2_.lock();
}

void SplitMergeEvent::setVertex2(VertexSPtr vertex_2) {
    this->vertex_2_ = vertex_2;
    this->neighborhood_2_ = VertexFacetNeighborhood(vertex_2);
}

FacetSPtr SplitMergeEvent::getFacet1() const {
    DEBUG_WPTR(facet_1_);
    return facet_1_.lock();
}

void SplitMergeEvent::setFacet1(FacetSPtr facet_1) {
    this->facet_1_ = facet_1;
}

FacetSPtr SplitMergeEvent::getFacet2() const {
    DEBUG_WPTR(facet_2_);
    return facet_2_.lock();
}

void SplitMergeEvent::setFacet2(FacetSPtr facet_2) {
    this->facet_2_ = facet_2;
}

void SplitMergeEvent::setHighlight(bool highlight) {
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

std::string SplitMergeEvent::toString() const {
    VertexSPtr vertex_1 = getVertex1();
    VertexSPtr vertex_2 = getVertex2();
    FacetSPtr facet_1 = getFacet1();
    FacetSPtr facet_2 = getFacet2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "SplitMergeEvent\n";
    sstr << "\t(ID=" << getID() << "; step ID=" << getStepID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(vertex1=" << vertex_1->toString() << ")\n";
    sstr << "\t(vertex2=" << vertex_2->toString() << ")\n";
    sstr << "\t(facet1=" << facet_1->getID() << ")\n";
    sstr << "\t(facet2=" << facet_2->getID() << ")";
    return sstr.str();
}

bool SplitMergeEvent::isValid() const {
    return node_ && !vertex_1_.expired() && !vertex_2_.expired() && !facet_1_.expired() && !facet_2_.expired();
}

bool SplitMergeEvent::isObsolete() const {
  if (VertexSPtr vertex_1 = getVertex1()) {
      // std::cout << "isObsolete(v" << vertex_1->getID() << ")?" << std::endl;
      if (!neighborhood_1_.checkNeighborhoodConsistency(vertex_1)) {
          return true;
      }
  }

  if (VertexSPtr vertex_2 = getVertex2()) {
      // std::cout << "isObsolete(v" << vertex_2->getID() << ")?" << std::endl;
      if (!neighborhood_2_.checkNeighborhoodConsistency(vertex_2)) {
          return true;
      }
  }

  return false;
}

bool SplitMergeEvent::operator==(const SplitMergeEvent& other) const {
    return (node_->getOffset() == other.node_->getOffset()) &&
           (*(node_->getPoint()) == *(other.node_->getPoint())) &&
           ((facet_1_.lock() == other.facet_1_.lock() &&
             facet_2_.lock() == other.facet_2_.lock()) ||
            (facet_1_.lock() == other.facet_2_.lock() &&
             facet_2_.lock() == other.facet_1_.lock())) &&
           ((vertex_1_.lock() == other.vertex_1_.lock() &&
             vertex_2_.lock() == other.vertex_2_.lock()) ||
            (vertex_1_.lock() == other.vertex_2_.lock() &&
             vertex_2_.lock() == other.vertex_1_.lock()));
}

} } }
