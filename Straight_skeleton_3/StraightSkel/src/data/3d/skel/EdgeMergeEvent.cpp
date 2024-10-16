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
 * @file   data/3d/skel/EdgeMergeEvent.cpp
 * @author Gernot Walzl
 * @date   2012-09-14
 */

#include "data/3d/skel/EdgeMergeEvent.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SkelEdgeData.h"
#include "util/StringFactory.h"

#include <sstream>

namespace data { namespace _3d { namespace skel {

EdgeMergeEvent::EdgeMergeEvent(PolyhedronSPtr polyhedron) : AbstractEvent(polyhedron) {
    type_ = AbstractEvent::EDGE_MERGE_EVENT;
}

EdgeMergeEvent::~EdgeMergeEvent() {
    node_.reset();
    facet_.reset();
    edge1_.reset();
    edge2_.reset();
}

EdgeMergeEventSPtr EdgeMergeEvent::create(PolyhedronSPtr polyhedron) {
    EdgeMergeEventSPtr result = EdgeMergeEventSPtr(new EdgeMergeEvent(polyhedron));
    return result;
}

NodeSPtr EdgeMergeEvent::getNode() const {
    DEBUG_SPTR(node_);
    return node_;
}

void EdgeMergeEvent::setNode(NodeSPtr node) {
    this->node_ = node;
}

CGAL::FT EdgeMergeEvent::getOffset() const {
    CGAL::FT result = 0.0;
    if (node_) {
        result = node_->getOffset();
    }
    return result;
}

FacetSPtr EdgeMergeEvent::getFacet() const {
    DEBUG_SPTR(facet_);
    return facet_;
}

void EdgeMergeEvent::setFacet(FacetSPtr facet) {
    this->facet_ = facet;
}

EdgeSPtr EdgeMergeEvent::getEdge1() const {
    DEBUG_SPTR(edge1_);
    return edge1_;
}

void EdgeMergeEvent::setEdge1(EdgeSPtr edge1) {
    this->edge1_ = edge1;
}

EdgeSPtr EdgeMergeEvent::getEdge2() const {
    DEBUG_SPTR(edge2_);
    return edge2_;
}

void EdgeMergeEvent::setEdge2(EdgeSPtr edge2) {
    this->edge2_ = edge2;
}

void EdgeMergeEvent::setHighlight(bool highlight) {
    if (!edge1_->hasData()) {
        SkelEdgeData::create(edge1_);
    }
    edge1_->getData()->setHighlight(highlight);
    if (!edge2_->hasData()) {
        SkelEdgeData::create(edge2_);
    }
    edge2_->getData()->setHighlight(highlight);
    EdgeSPtr edge_toremove_1 = edge1_->next(facet_);
    EdgeSPtr edge_toremove_2 = edge_toremove_1->next(facet_);
    if (!edge_toremove_1->hasData()) {
        SkelEdgeData::create(edge_toremove_1);
    }
    edge_toremove_1->getData()->setHighlight(highlight);
    if (!edge_toremove_2->hasData()) {
        SkelEdgeData::create(edge_toremove_2);
    }
    edge_toremove_2->getData()->setHighlight(highlight);
}

std::string EdgeMergeEvent::toString() const {
    std::stringstream sstr;
    sstr.precision(17);
    sstr << "EdgeMergeEvent\n";
    sstr << "\t(offset=" << util::StringFactory::fromDouble(CGAL::to_double(getOffset())) << ")\n";
    sstr << "\t(node=" << *(getNode()->getPoint()) << ")\n";
    sstr << "\t(facet=" << getFacet()->getID() << "\n";
    sstr << "\t(edgeA=" << getEdge1()->getID() << "\n\t\t[" << getEdge1()->getVertexSrc()->toString() << "\n\t\t "
                                                            << getEdge1()->getVertexDst()->toString() << "])\n";
    sstr << "\t(edgeB=" << getEdge2()->getID() << "\n\t\t[" << getEdge2()->getVertexSrc()->toString() << "\n\t\t "
                                                            << getEdge2()->getVertexDst()->toString() << "])\n";
    return sstr.str();
}

} } }
