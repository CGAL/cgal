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
 * @file   data/3d/skel/SurfaceEvent.cpp
 * @author Gernot Walzl
 * @date   2012-09-10
 */

#include "data/3d/skel/SurfaceEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "util/StringFactory.h"

#include <sstream>
#include <string>

namespace data { namespace _3d { namespace skel {

SurfaceEvent::SurfaceEvent() {
    type_ = AbstractEvent::SURFACE_EVENT;
}

SurfaceEvent::~SurfaceEvent() {
    node_.reset();
    edge1_.reset();
    edge2_.reset();
}

SurfaceEventSPtr SurfaceEvent::create() {
    SurfaceEventSPtr result = SurfaceEventSPtr(new SurfaceEvent());
    return result;
}

NodeSPtr SurfaceEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void SurfaceEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT SurfaceEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

EdgeSPtr SurfaceEvent::getEdge1() const {
    DEBUG_WPTR(edge1_);
    return edge1_.lock();
}

void SurfaceEvent::setEdge1(EdgeSPtr edge1) {
    this->edge1_ = edge1;
}

EdgeSPtr SurfaceEvent::getEdge2() const {
    DEBUG_WPTR(edge2_);
    return edge2_.lock();
}

void SurfaceEvent::setEdge2(EdgeSPtr edge2) {
    this->edge2_ = edge2;
}

void SurfaceEvent::setHighlight(bool highlight) {
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

std::string SurfaceEvent::toString() const {
    EdgeSPtr edge1 = getEdge1();
    EdgeSPtr edge2 = getEdge2();

    std::stringstream sstr;
    sstr.precision(17);
    sstr << "SurfaceEvent\n";
    sstr << "\t(ID=" << getID() << ")\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(edgeA=" << edge1->getID() << "\n\t\t[" << edge1->getVertexSrc()->toString() << "\n\t\t "
                                                       << edge1->getVertexDst()->toString() << "]\n"
         << "\t edgeB=" << edge2->getID() << "\n\t\t[" << edge2->getVertexSrc()->toString() << "\n\t\t "
                                                       << edge2->getVertexDst()->toString() << "])";

    return sstr.str();
}

bool SurfaceEvent::isValid() const {
    return node_ && !edge1_.expired() && !edge2_.expired();
}

} } }
