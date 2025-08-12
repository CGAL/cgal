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
 * @file   data/3d/skel/SkelVertexData.cpp
 * @author Gernot Walzl
 * @date   2012-04-05
 */

#include "data/3d/skel/SkelVertexData.h"

#include "debug.h"
#include "data/3d/Vertex.h"

namespace data { namespace _3d { namespace skel {

SkelVertexData::SkelVertexData() {
    // intentionally does nothing
}

SkelVertexData::~SkelVertexData() {
    // intentionally does nothing
}

SkelVertexDataSPtr SkelVertexData::create(VertexSPtr vertex) {
    SkelVertexDataSPtr result = SkelVertexDataSPtr(new SkelVertexData());
    result->setVertex(vertex);
    vertex->setData(result);
    return result;
}

ArcSPtr SkelVertexData::getArc() const {
    CGAL_SS3_DEBUG_WPTR(arc_);
    return this->arc_.lock();
}

void SkelVertexData::setArc(ArcSPtr arc) {
    this->arc_ = arc;
}

NodeSPtr SkelVertexData::getNode() const {
    return this->node_.lock();
}

void SkelVertexData::setNode(NodeSPtr node) {
    this->node_ = node;
}

VertexSPtr SkelVertexData::getOffsetVertex() const {
    CGAL_SS3_DEBUG_WPTR(offset_vertex_);
    return this->offset_vertex_.lock();
}

void SkelVertexData::setOffsetVertex(VertexSPtr offset_vertex) {
    this->offset_vertex_ = offset_vertex;
}

} } }
