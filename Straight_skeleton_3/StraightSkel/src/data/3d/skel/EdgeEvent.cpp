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
 * @file   data/3d/skel/EdgeEvent.cpp
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#include "data/3d/skel/EdgeEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "util/StringFactory.h"

#include <sstream>
#include <string>

namespace data { namespace _3d { namespace skel {

EdgeEvent::EdgeEvent() {
    type_ = AbstractEvent::EDGE_EVENT;
}

EdgeEvent::~EdgeEvent() {
    node_.reset();
    edge_.reset(); // @fixme is there still a point since edge_ is a weak pointer?
}

EdgeEventSPtr EdgeEvent::create() {
    EdgeEventSPtr result = EdgeEventSPtr(new EdgeEvent());
    return result;
}

NodeSPtr EdgeEvent::getNode() const {
    CGAL_SS3_DEBUG_SPTR(node_);
    return node_;
}

void EdgeEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT EdgeEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr EdgeEvent::getEdge() const {
    CGAL_SS3_DEBUG_WPTR(edge_);
    return edge_.lock();
}

void EdgeEvent::setEdge(EdgeSPtr edge) {
    this->edge_ = edge;
    this->neighborhood_ = EdgeFacetNeighborhood(edge);
}

void EdgeEvent::setHighlight(bool highlight) {
    EdgeSPtr edge = getEdge();
    if (!edge->hasData()) {
        SkelEdgeData::create(edge);
    }
    edge->getData()->setHighlight(highlight);
}

std::string EdgeEvent::toString() const {
    EdgeSPtr edge = getEdge();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "EdgeEvent\n";
    sstr << "\t(ID=" << getID() << "; step ID=" << getStepID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(edgeA=" << edge->getID() << "\n\t\t[" << edge->getVertexSrc()->toString() << "\n\t\t "
                                                      << edge->getVertexDst()->toString() << "])";
    return sstr.str();
}

bool EdgeEvent::isValid() const {
    return node_ && !edge_.expired();
}

bool EdgeEvent::isObsolete() const {
    if (EdgeSPtr edge = getEdge()) {
        // std::cout << "isObsolete(e " << edge->getID() << ")?" << std::endl;
        return ! neighborhood_.checkNeighborhoodConsistency(edge);
    }

    return false;
}

bool EdgeEvent::operator==(const EdgeEvent& other) const {
    return (node_->getOffset() == other.node_->getOffset()) &&
           (*(node_->getPoint()) == *(other.node_->getPoint())) &&
           (edge_.lock() == other.edge_.lock());
}


} } }
