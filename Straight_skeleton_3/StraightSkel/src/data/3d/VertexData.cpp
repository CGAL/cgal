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
 * @file   data/3d/VertexData.cpp
 * @author Gernot Walzl
 * @date   2011-11-22
 */

#include "data/3d/VertexData.h"

#include "data/3d/Vertex.h"

namespace data { namespace _3d {

VertexData::VertexData() {
    highlight_ = false;
}

VertexData::~VertexData() {
    // intentionally does nothing
}

VertexDataSPtr VertexData::create(VertexSPtr vertex) {
    VertexDataSPtr result = VertexDataSPtr(new VertexData());
    result->setVertex(vertex);
    vertex->setData(result);
    return result;
}

VertexSPtr VertexData::getVertex() const {
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
