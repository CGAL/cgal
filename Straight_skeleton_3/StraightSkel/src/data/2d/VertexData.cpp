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
    DEBUG_WPTR(vertex_);
    if (this->vertex_.expired())
        return VertexSPtr();
    else
        return VertexSPtr(this->vertex_);
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
