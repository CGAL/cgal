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
 * @file   data/3d/CircularVertexData.cpp
 * @author Gernot Walzl
 * @date   2012-11-30
 */

#include "data/3d/CircularVertexData.h"

#include "debug.h"
#include "data/3d/CircularVertex.h"

namespace data { namespace _3d {

CircularVertexData::CircularVertexData() {
    highlight_ = false;
}

CircularVertexData::~CircularVertexData() {
    // intentionally does nothing
}

CircularVertexSPtr CircularVertexData::getVertex() const {
    CGAL_SS3_DEBUG_WPTR(vertex_);
    return this->vertex_.lock();
}

void CircularVertexData::setVertex(CircularVertexSPtr vertex) {
    this->vertex_ = vertex;
}

bool CircularVertexData::isHighlight() const {
    return highlight_;
}

void CircularVertexData::setHighlight(bool highlight) {
    highlight_ = highlight;
}

} }
