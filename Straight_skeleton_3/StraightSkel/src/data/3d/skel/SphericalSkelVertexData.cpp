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
    //DEBUG_WPTR(arc_);
    if (this->arc_.expired())
        return CircularArcSPtr();
    else
        return CircularArcSPtr(this->arc_);
}

void SphericalSkelVertexData::setArc(CircularArcSPtr arc) {
    this->arc_ = arc;
}

CircularNodeSPtr SphericalSkelVertexData::getNode() const {
    DEBUG_WPTR(node_);
    if (this->node_.expired())
        return CircularNodeSPtr();
    else
        return CircularNodeSPtr(this->node_);
}

void SphericalSkelVertexData::setNode(CircularNodeSPtr node) {
    this->node_ = node;
}

CircularVertexSPtr SphericalSkelVertexData::getOffsetVertex() const {
    DEBUG_WPTR(offset_vertex_);
    if (this->offset_vertex_.expired())
        return CircularVertexSPtr();
    else
        return CircularVertexSPtr(this->offset_vertex_);
}

void SphericalSkelVertexData::setOffsetVertex(CircularVertexSPtr offset_vertex) {
    this->offset_vertex_ = offset_vertex;
}

double SphericalSkelVertexData::getSpeed() const {
    return this->speed_;
}

void SphericalSkelVertexData::setSpeed(double speed) {
    this->speed_ = speed;
}

EdgeSPtr SphericalSkelVertexData::getEdgeOrigin() const {
    DEBUG_WPTR(edge_origin_);
    if (this->edge_origin_.expired())
        return EdgeSPtr();
    else
        return EdgeSPtr(this->edge_origin_);
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
