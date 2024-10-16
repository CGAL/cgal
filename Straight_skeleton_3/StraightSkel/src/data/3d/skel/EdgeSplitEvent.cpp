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
 * @file   data/3d/skel/EdgeSplitEvent.cpp
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#include "data/3d/skel/EdgeSplitEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "util/StringFactory.h"

#include <sstream>
#include <string>

namespace data { namespace _3d { namespace skel {

EdgeSplitEvent::EdgeSplitEvent(PolyhedronSPtr polyhedron) : AbstractEvent(polyhedron) {
    type_ = AbstractEvent::EDGE_SPLIT_EVENT;
}

EdgeSplitEvent::~EdgeSplitEvent() {
    node_.reset();
    edge1_.reset();
    edge2_.reset();
}

EdgeSplitEventSPtr EdgeSplitEvent::create(PolyhedronSPtr polyhedron) {
    EdgeSplitEventSPtr result = EdgeSplitEventSPtr(new EdgeSplitEvent(polyhedron));
    return result;
}

NodeSPtr EdgeSplitEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void EdgeSplitEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT EdgeSplitEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr EdgeSplitEvent::getEdge1() const {
    DEBUG_SPTR(edge1_);
    return edge1_;
}

void EdgeSplitEvent::setEdge1(EdgeSPtr edge1) {
    this->edge1_ = edge1;
}

EdgeSPtr EdgeSplitEvent::getEdge2() const {
    DEBUG_SPTR(edge2_);
    return edge2_;
}

void EdgeSplitEvent::setEdge2(EdgeSPtr edge2) {
    this->edge2_ = edge2;
}

void EdgeSplitEvent::setHighlight(bool highlight) {
    if (!edge1_->hasData()) {
        SkelEdgeData::create(edge1_);
    }
    edge1_->getData()->setHighlight(highlight);
    if (!edge2_->hasData()) {
        SkelEdgeData::create(edge2_);
    }
    edge2_->getData()->setHighlight(highlight);
}

std::string EdgeSplitEvent::toString() const {
    std::stringstream sstr;
    sstr.precision(17);
    sstr << "EdgeSplitEvent\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(edgeA=" << getEdge1()->getID() << "[" << getEdge1()->getVertexSrc()->getID() << "-"
                                                      << getEdge1()->getVertexDst()->getID() << "]"
         << "; edgeB=" << getEdge2()->getID() << "[" << getEdge2()->getVertexSrc()->getID() << "-"
                                                     << getEdge2()->getVertexDst()->getID() << "]";
    return sstr.str();
}

} } }
