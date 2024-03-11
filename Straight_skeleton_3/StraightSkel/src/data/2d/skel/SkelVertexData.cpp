/**
 * @file   data/2d/skel/SkelVertexData.cpp
 * @author Gernot Walzl
 * @date   2012-02-09
 */

#include "data/2d/skel/SkelVertexData.h"

#include "data/2d/Vertex.h"

namespace data { namespace _2d { namespace skel {

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
    DEBUG_WPTR(arc_);
    if (this->arc_.expired())
        return ArcSPtr();
    else
        return ArcSPtr(this->arc_);
}

void SkelVertexData::setArc(ArcSPtr arc) {
    this->arc_ = arc;
}

NodeSPtr SkelVertexData::getNode() const {
    DEBUG_WPTR(node_);
    if (this->node_.expired())
        return NodeSPtr();
    else
        return NodeSPtr(this->node_);
}

void SkelVertexData::setNode(NodeSPtr node) {
    this->node_ = node;
}

VertexSPtr SkelVertexData::getOffsetVertex() const {
    DEBUG_WPTR(offset_vertex_);
    if (this->offset_vertex_.expired())
        return VertexSPtr();
    else
        return VertexSPtr(this->offset_vertex_);
}

void SkelVertexData::setOffsetVertex(VertexSPtr offset_vertex) {
    this->offset_vertex_ = offset_vertex;
}

} } }
