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
 * @file   data/2d/skel/SkelEdgeData.cpp
 * @author Gernot Walzl
 * @date   2012-05-12
 */

#include "data/2d/skel/SkelEdgeData.h"

#include "data/2d/Edge.h"
#include "debug.h"

namespace data { namespace _2d { namespace skel {

SkelEdgeData::SkelEdgeData() {
    speed_ = 1.0;
}

SkelEdgeData::~SkelEdgeData() {
    // intentionally does nothing
}

SkelEdgeDataSPtr SkelEdgeData::create(EdgeSPtr edge) {
    SkelEdgeDataSPtr result = SkelEdgeDataSPtr(new SkelEdgeData());
    result->setEdge(edge);
    result->setEdgeOrigin(edge);
    edge->setData(result);
    return result;
}

EdgeSPtr SkelEdgeData::getOffsetEdge() const {
    CGAL_SS3_DEBUG_WPTR(offset_edge_);
    if (this->offset_edge_.expired())
        return EdgeSPtr();
    else
        return EdgeSPtr(this->offset_edge_);
}

void SkelEdgeData::setOffsetEdge(EdgeSPtr offset_edge) {
    this->offset_edge_ = offset_edge;
}

EdgeSPtr SkelEdgeData::getEdgeOrigin() const {
    CGAL_SS3_DEBUG_WPTR(edge_origin_);
    if (this->edge_origin_.expired())
        return EdgeSPtr();
    else
        return EdgeSPtr(this->edge_origin_);
}

void SkelEdgeData::setEdgeOrigin(EdgeSPtr edge_origin) {
    this->edge_origin_ = edge_origin;
}

CGAL::FT SkelEdgeData::getSpeed() const {
    return speed_;
}

void SkelEdgeData::setSpeed(CGAL::FT speed) {
    speed_ = speed;
}

} } }
