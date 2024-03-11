/**
 * @file   data/3d/CircularEdge.cpp
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#include "data/3d/CircularEdge.h"

#include "debug.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/SphericalPolygon.h"
#include "util/StringFactory.h"

namespace data { namespace _3d {

CircularEdge::CircularEdge(CircularVertexSPtr src, CircularVertexSPtr dst) {
    vertex_src_ = src;
    vertex_dst_ = dst;
}

CircularEdge::~CircularEdge() {
    vertex_src_.reset();
    vertex_dst_.reset();
}

CircularEdgeSPtr CircularEdge::create(CircularVertexSPtr src, CircularVertexSPtr dst) {
    CircularEdgeSPtr result = CircularEdgeSPtr(new CircularEdge(src, dst));
    src->setEdgeOut(result);
    dst->setEdgeIn(result);
    return result;
}

CircularVertexSPtr CircularEdge::getVertexSrc() const {
    DEBUG_SPTR(vertex_src_);
    return this->vertex_src_;
}

void CircularEdge::setVertexSrc(CircularVertexSPtr src) {
    this->vertex_src_ = src;
}

CircularVertexSPtr CircularEdge::getVertexDst() const {
    DEBUG_SPTR(vertex_dst_);
    return this->vertex_dst_;
}

void CircularEdge::setVertexDst(CircularVertexSPtr dst) {
    this->vertex_dst_ = dst;
}

SphericalPolygonSPtr CircularEdge::getPolygon() const {
    DEBUG_WPTR(polygon_);
    if (this->polygon_.expired())
        return SphericalPolygonSPtr();
    else
        return SphericalPolygonSPtr(this->polygon_);
}

void CircularEdge::setPolygon(SphericalPolygonSPtr polygon) {
    this->polygon_ = polygon;
}

std::list<CircularEdgeSPtr>::iterator CircularEdge::getListIt() const {
    return this->list_it_;
}

void CircularEdge::setListIt(std::list<CircularEdgeSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

CircularEdgeDataSPtr CircularEdge::getData() const {
    DEBUG_SPTR(data_);
    return this->data_;
}

void CircularEdge::setData(CircularEdgeDataSPtr data) {
    this->data_ = data;
}

bool CircularEdge::hasData() const {
    bool result = false;
    if (data_) {
        result = true;
    }
    return result;
}

Plane3SPtr CircularEdge::getSupportingPlane() const {
    DEBUG_SPTR(supporting_plane_);
    return supporting_plane_;
}

void CircularEdge::setSupportingPlane(Plane3SPtr supporting_plane) {
    this->supporting_plane_ = supporting_plane;
}

bool CircularEdge::initSupportingPlane() {
    bool result = false;
    if (!polygon_.expired()) {
        Sphere3SPtr sphere = getPolygon()->getSphere();
        Point3SPtr p_center = KernelFactory::createPoint3(sphere);
        Point3SPtr p_src = vertex_src_->getPoint();
        Point3SPtr p_dst = vertex_dst_->getPoint();
        supporting_plane_ = KernelFactory::createPlane3(p_center, p_src, p_dst);
        if (supporting_plane_) {
            result = true;
        }
    }
    return result;
}

Plane3SPtr CircularEdge::supportingPlane() {
    if (!supporting_plane_) {
        initSupportingPlane();
    }
    DEBUG_SPTR(supporting_plane_);
    return supporting_plane_;
}

CircularEdgeSPtr CircularEdge::next() const {
    return vertex_dst_->getEdgeOut();
}

CircularEdgeSPtr CircularEdge::prev() const {
    return vertex_src_->getEdgeIn();
}

std::string CircularEdge::toString() const {
    std::string result("CircularEdge(");
    result += util::StringFactory::fromPointer(this) + ", ";
    result += "src=" + vertex_src_->toString() + ", ";
    result += "dst=" + vertex_dst_->toString() + ")";
    return result;
}

} }
