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

EdgeSplitEvent::EdgeSplitEvent() {
    type_ = AbstractEvent::EDGE_SPLIT_EVENT;
}

EdgeSplitEvent::~EdgeSplitEvent() {
    node_.reset();
    edge1_.reset();
    edge2_.reset();
}

EdgeSplitEventSPtr EdgeSplitEvent::create() {
    EdgeSplitEventSPtr result = EdgeSplitEventSPtr(new EdgeSplitEvent());
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
    DEBUG_WPTR(edge1_);
    return edge1_.lock();
}

void EdgeSplitEvent::setEdge1(EdgeSPtr edge1) {
    this->edge1_ = edge1;
    this->neighborhood1_ = EdgeFacetNeighborhood(edge1);
}

EdgeSPtr EdgeSplitEvent::getEdge2() const {
    DEBUG_WPTR(edge2_);
    return edge2_.lock();
}

void EdgeSplitEvent::setEdge2(EdgeSPtr edge2) {
    this->edge2_ = edge2;
    this->neighborhood2_ = EdgeFacetNeighborhood(edge2);
}

int EdgeSplitEvent::getEdgeOrientation() const {
    return edge_orientation_;
}

void EdgeSplitEvent::setEdgeOrientation(int edge_orientation) {
    this->edge_orientation_ = edge_orientation;
}

void EdgeSplitEvent::setHighlight(bool highlight) {
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

std::string EdgeSplitEvent::toString() const {
    EdgeSPtr edge1 = getEdge1();
    EdgeSPtr edge2 = getEdge2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "EdgeSplitEvent\n";
    sstr << "\t(ID=" << getID() << "; step ID=" << getStepID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(edgeA=" << edge1->getID() << "[" << edge1->getVertexSrc()->getID() << "-"
                                                 << edge1->getVertexDst()->getID() << "]"
         << "; edgeB=" << edge2->getID() << "[" << edge2->getVertexSrc()->getID() << "-"
                                                << edge2->getVertexDst()->getID() << "])\n";
    sstr << "\t(edge_orientation=" << edge_orientation_ << ")\n";
    return sstr.str();
}

bool EdgeSplitEvent::isValid() const {
    return node_ && !edge1_.expired() && !edge2_.expired();
}

bool EdgeSplitEvent::isObsolete() const {
    if (EdgeSPtr edge_1 = getEdge1()) {
        // std::cout << "isObsolete(e " << edge_1->getID() << ")?" << std::endl;
        if (!neighborhood1_.checkNeighborhoodConsistency(edge_1)) {
            return true;
        }
    }

    if (EdgeSPtr edge_2 = getEdge2()) {
        // std::cout << "isObsolete(e " << edge_2->getID() << ")?" << std::endl;
        if (!neighborhood2_.checkNeighborhoodConsistency(edge_2)) {
            return true;
        }
    }

    return false;
}

bool EdgeSplitEvent::operator==(const EdgeSplitEvent& other) const {
    return (node_->getOffset() == other.node_->getOffset()) &&
           (*(node_->getPoint()) == *(other.node_->getPoint())) &&
           ((edge1_.lock() == other.edge1_.lock() &&
             edge2_.lock() == other.edge2_.lock()) ||
            (edge1_.lock() == other.edge2_.lock() &&
             edge2_.lock() == other.edge1_.lock())) &&
           (edge_orientation_ == other.edge_orientation_);
}

} } }
