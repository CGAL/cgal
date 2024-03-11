/**
 * @file   data/2d/Edge.cpp
 * @author Gernot Walzl
 * @date   2011-11-22
 */

#include "data/2d/Edge.h"

#include "debug.h"
#include "util/StringFactory.h"
#include "data/2d/Vertex.h"
#include "data/2d/KernelFactory.h"

namespace data { namespace _2d {

Edge::Edge() {
    this->id_ = -1;
}

Edge::~Edge() {
    // intentionally does nothing
}

EdgeSPtr Edge::create(VertexSPtr src, VertexSPtr dst) {
    EdgeSPtr result = EdgeSPtr(new Edge());
    result->setVertexSrc(src);
    result->setVertexDst(dst);
    src->setEdgeOut(result);
    dst->setEdgeIn(result);
    return result;
}

VertexSPtr Edge::getVertexSrc() const {
    DEBUG_SPTR(vertex_src_);
    return this->vertex_src_;
}

void Edge::setVertexSrc(VertexSPtr src) {
    this->vertex_src_ = src;
}

VertexSPtr Edge::getVertexDst() const {
    DEBUG_SPTR(vertex_dst_);
    return this->vertex_dst_;
}

void Edge::setVertexDst(VertexSPtr dst) {
    this->vertex_dst_ = dst;
}

PolygonSPtr Edge::getPolygon() const {
    DEBUG_WPTR(polygon_);
    if (this->polygon_.expired())
        return PolygonSPtr();
    else
        return PolygonSPtr(this->polygon_);
}

void Edge::setPolygon(PolygonSPtr polygon) {
    this->polygon_ = polygon;
}

std::list<EdgeSPtr>::iterator Edge::getListIt() const {
    return this->list_it_;
}

void Edge::setListIt(std::list<EdgeSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

EdgeDataSPtr Edge::getData() const {
    DEBUG_SPTR(data_);
    return this->data_;
}

void Edge::setData(EdgeDataSPtr data) {
    this->data_ = data;
}

bool Edge::hasData() const {
    bool result = false;
    if (data_) {
        result = true;
    }
    return result;
}

int Edge::getID() const {
    return this->id_;
}

void Edge::setID(int id) {
    this->id_ = id;
}

EdgeSPtr Edge::next() const {
    return vertex_dst_->getEdgeOut();
}

EdgeSPtr Edge::prev() const {
    return vertex_src_->getEdgeIn();
}

Segment2SPtr Edge::segment() const {
    return KernelFactory::createSegment2(
            vertex_src_->getPoint(), vertex_dst_->getPoint());
}

Line2SPtr Edge::line() const {
    return KernelFactory::createLine2(
            vertex_src_->getPoint(), vertex_dst_->getPoint());
}

std::string Edge::toString() const {
    std::string result("Edge(");
    if (id_ != -1) {
        result += "id=" + util::StringFactory::fromInteger(id_) + ", ";
    } else {
        result += util::StringFactory::fromPointer(this) + ", ";
    }
    result += "src=" + vertex_src_->toString() + ", ";
    result += "dst=" + vertex_dst_->toString() + ")";
    return result;
}

} }
