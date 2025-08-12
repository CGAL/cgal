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
 * @file   data/3d/skel/SphericalSkelVertexData.cpp
 * @author Gernot Walzl
 * @date   2012-11-30
 */

#include "data/3d/skel/SphericalSkelVertexData.h"

#include "debug.h"
#include "data/3d/CircularVertex.h"

namespace data { namespace _3d { namespace skel {

SphericalSkelVertexData::SphericalSkelVertexData() {
    this->speed_ = 1.0;
    this->reflex_ = false;
    this->invert_ = true;
}

SphericalSkelVertexData::~SphericalSkelVertexData() {
    // intentionally does nothing
}

SphericalSkelVertexDataSPtr SphericalSkelVertexData::create(CircularVertexSPtr vertex) {
    SphericalSkelVertexDataSPtr result = SphericalSkelVertexDataSPtr(new SphericalSkelVertexData());
    result->setVertex(vertex);
    vertex->setData(result);
    return result;
}

CircularArcSPtr SphericalSkelVertexData::getArc() const {
    //CGAL_SS3_DEBUG_WPTR(arc_);
    return this->arc_.lock();
}

void SphericalSkelVertexData::setArc(CircularArcSPtr arc) {
    this->arc_ = arc;
}

CircularNodeSPtr SphericalSkelVertexData::getNode() const {
    CGAL_SS3_DEBUG_WPTR(node_);
    return this->node_.lock();
}

void SphericalSkelVertexData::setNode(CircularNodeSPtr node) {
    this->node_ = node;
}

CircularVertexSPtr SphericalSkelVertexData::getOffsetVertex() const {
    CGAL_SS3_DEBUG_WPTR(offset_vertex_);
    return this->offset_vertex_.lock();
}

void SphericalSkelVertexData::setOffsetVertex(CircularVertexSPtr offset_vertex) {
    this->offset_vertex_ = offset_vertex;
}

CGAL::FT SphericalSkelVertexData::getSpeed() const {
    return this->speed_;
}

void SphericalSkelVertexData::setSpeed(CGAL::FT speed) {
    this->speed_ = speed;
}

EdgeSPtr SphericalSkelVertexData::getEdgeOrigin() const {
    CGAL_SS3_DEBUG_WPTR(edge_origin_);
    return this->edge_origin_.lock();
}

void SphericalSkelVertexData::setEdgeOrigin(EdgeSPtr edge_origin) {
    this->edge_origin_ = edge_origin;
}

bool SphericalSkelVertexData::isReflex() const {
    return this->reflex_;
}

void SphericalSkelVertexData::setReflex(bool reflex) {
    this->reflex_ = reflex;
}

bool SphericalSkelVertexData::doInvert() const {
    return this->invert_;
}

void SphericalSkelVertexData::setInvert(bool invert) {
    this->invert_ = invert;
}

} } }
