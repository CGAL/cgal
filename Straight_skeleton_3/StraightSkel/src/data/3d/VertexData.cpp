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
