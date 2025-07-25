// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#include "data/3d/skel/VanishEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "util/StringFactory.h"

#include <sstream>
#include <string>

namespace data { namespace _3d { namespace skel {

VanishEvent::VanishEvent() {
    type_ = AbstractEvent::VANISH_EVENT;
}

VanishEvent::~VanishEvent() {
    node_.reset();
    edge_.reset(); // @fixme is there still a point since edge_ is a weak pointer?
}

VanishEventSPtr VanishEvent::create() {
    VanishEventSPtr result = VanishEventSPtr(new VanishEvent());
    return result;
}

NodeSPtr VanishEvent::getNode() const {
    CGAL_SS3_DEBUG_SPTR(node_);
    return node_;
}

void VanishEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT VanishEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr VanishEvent::getEdge() const {
    CGAL_SS3_DEBUG_WPTR(edge_);
    return edge_.lock();
}

void VanishEvent::setEdge(EdgeSPtr edge) {
    this->edge_ = edge;
    this->neighborhood_ = EdgeFacetNeighborhood(edge);
}

void VanishEvent::setHighlight(bool highlight) {
    EdgeSPtr edge = getEdge();
    if (!edge->hasData()) {
        SkelEdgeData::create(edge);
    }
    edge->getData()->setHighlight(highlight);
}

std::string VanishEvent::toString() const {
    EdgeSPtr edge = getEdge();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "VanishEvent\n";
    sstr << "\t(ID=" << getID() << "; step ID=" << getStepID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(edgeA=" << edge->getID() << "\n\t\t[" << edge->getVertexSrc()->toString() << "\n\t\t "
                                                      << edge->getVertexDst()->toString() << "])";
    return sstr.str();
}

bool VanishEvent::isValid() const {
    return node_ && !edge_.expired();
}


bool VanishEvent::isObsolete() const {
    if (EdgeSPtr edge = getEdge()) {
        return ! neighborhood_.checkNeighborhoodConsistency(edge);
    }

    return false;
}

bool VanishEvent::operator==(const VanishEvent& other) const {
    return (node_->getOffset() == other.node_->getOffset()) &&
           // && (edge_.lock() == other.edge_.lock()) // because of multiple reps...
           (*(node_->getPoint()) == *(other.node_->getPoint()));
}

} } }
