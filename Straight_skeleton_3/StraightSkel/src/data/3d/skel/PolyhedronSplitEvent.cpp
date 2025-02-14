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
 * @file   data/3d/skel/PolyhedronSplitEvent.cpp
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#include "data/3d/skel/PolyhedronSplitEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "util/StringFactory.h"

#include <sstream>
#include <string>

namespace data { namespace _3d { namespace skel {

PolyhedronSplitEvent::PolyhedronSplitEvent() {
    type_ = AbstractEvent::POLYHEDRON_SPLIT_EVENT;
}

PolyhedronSplitEvent::~PolyhedronSplitEvent() {
    node_.reset();
    edge1_.reset();
    edge2_.reset();
}

PolyhedronSplitEventSPtr PolyhedronSplitEvent::create() {
    PolyhedronSplitEventSPtr result = PolyhedronSplitEventSPtr(new PolyhedronSplitEvent());
    return result;
}

NodeSPtr PolyhedronSplitEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void PolyhedronSplitEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT PolyhedronSplitEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr PolyhedronSplitEvent::getEdge1() const {
    DEBUG_WPTR(edge1_);
    return edge1_.lock();
}

void PolyhedronSplitEvent::setEdge1(EdgeSPtr edge1) {
    this->edge1_ = edge1;
}

EdgeSPtr PolyhedronSplitEvent::getEdge2() const {
    DEBUG_WPTR(edge2_);
    return edge2_.lock();
}

void PolyhedronSplitEvent::setEdge2(EdgeSPtr edge2) {
    this->edge2_ = edge2;
}

void PolyhedronSplitEvent::setHighlight(bool highlight) {
    EdgeSPtr edge1 = getEdge1();
    EdgeSPtr edge2 = getEdge2();

    if (!edge1->hasData()) {
        SkelEdgeData::create(edge1);
    }
    edge1->getData()->setHighlight(highlight);
    if (!edge2->hasData()) {
        SkelEdgeData::create(edge2);
    }
    edge2->getData()->setHighlight(highlight);
}

std::string PolyhedronSplitEvent::toString() const {
    EdgeSPtr edge1 = getEdge1();
    EdgeSPtr edge2 = getEdge2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "PolyhedronSplitEvent\n";
    sstr << "\t(ID=" << getID() << "; step ID=" << getStepID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(edgeA=" << edge1->getID() << "[" << edge1->getVertexSrc()->getID() << "-"
                                                 << edge1->getVertexDst()->getID() << "]"
         << "; edgeB=" << edge2->getID() << "[" << edge2->getVertexSrc()->getID() << "-"
                                                << edge2->getVertexDst()->getID() << "])";
    return sstr.str();
}

bool PolyhedronSplitEvent::isValid() const {
    return node_ && !edge1_.expired() && !edge2_.expired();
}

} } }
