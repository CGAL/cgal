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
 * @file   data/2d/VertexData.cpp
 * @author Gernot Walzl
 * @date   2011-11-22
 */

#include "data/2d/VertexData.h"

#include "debug.h"

namespace data { namespace _2d {

VertexData::VertexData() {
    highlight_ = false;
}

VertexData::~VertexData() {
    // intentionally does nothing
}

VertexSPtr VertexData::getVertex() const {
    CGAL_SS3_DEBUG_WPTR(vertex_);
    return this->vertex_.lock();
}

void VertexData::setVertex(VertexSPtr vertex) {
    this->vertex_ = vertex;
}

bool VertexData::isHighlight() const {
    return highlight_;
}

void VertexData::setHighlight(bool highlight) {
    highlight_ = highlight;
}

} }
