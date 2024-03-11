/**
 * @file   data/3d/skel/CircularArc.cpp
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#include "data/3d/skel/CircularArc.h"

#include "debug.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/skel/CircularNode.h"
#include "data/3d/skel/SphericalSkeleton.h"
#include "util/StringFactory.h"

namespace data { namespace _3d { namespace skel {

CircularArc::CircularArc(CircularNodeSPtr node_src, Vector3SPtr direction) {
    node_src_ = node_src;
    node_dst_ = CircularNodeSPtr();
    direction_ = direction;
}

CircularArc::CircularArc(CircularNodeSPtr node_src, CircularNodeSPtr node_dst) {
    node_src_ = node_src;
    node_dst_ = node_dst;
    direction_ = Vector3SPtr();
}

CircularArc::~CircularArc() {
    node_src_.reset();
    node_dst_.reset();
}

CircularArcSPtr CircularArc::create(CircularNodeSPtr node_src, Vector3SPtr direction) {
    CircularArcSPtr result = CircularArcSPtr(new CircularArc(node_src, direction));
    node_src->addArc(result);
    return result;
}

CircularArcSPtr CircularArc::create(CircularNodeSPtr node_src, CircularNodeSPtr node_dst) {
    CircularArcSPtr result = CircularArcSPtr(new CircularArc(node_src, node_dst));
    node_src->addArc(result);
    node_dst->addArc(result);
    return result;
}

CircularNodeSPtr CircularArc::getNodeSrc() const {
    DEBUG_SPTR(node_src_);
    return node_src_;
}

void CircularArc::setNodeSrc(CircularNodeSPtr node_src) {
    this->node_src_ = node_src;
}


std::list<CircularArcWPtr>::iterator CircularArc::getNodeSrcListIt() const {
    return this->node_src_list_it_;
}

void CircularArc::setNodeSrcListIt(std::list<CircularArcWPtr>::iterator node_src_list_it) {
    this->node_src_list_it_ = node_src_list_it;
}

CircularNodeSPtr CircularArc::getNodeDst() const {
    DEBUG_SPTR(node_dst_);
    return node_dst_;
}

void CircularArc::setNodeDst(CircularNodeSPtr node_dst) {
    this->node_dst_ = node_dst;
}

std::list<CircularArcWPtr>::iterator CircularArc::getNodeDstListIt() const {
    return this->node_dst_list_it_;
}

void CircularArc::setNodeDstListIt(std::list<CircularArcWPtr>::iterator node_dst_list_it) {
    this->node_dst_list_it_ = node_dst_list_it;
}

Vector3SPtr CircularArc::getDirection() const {
    DEBUG_SPTR(direction_);
    return this->direction_;
}

void CircularArc::setDirection(Vector3SPtr direction) {
    this->direction_ = direction;
}

CircularEdgeSPtr CircularArc::getEdgeLeft() const {
    DEBUG_SPTR(edge_left_);
    return edge_left_;
}

void CircularArc::setEdgeLeft(CircularEdgeSPtr edge_left) {
    this->edge_left_ = edge_left;
}

CircularEdgeSPtr CircularArc::getEdgeRight() const {
    DEBUG_SPTR(edge_right_);
    return edge_right_;
}

void CircularArc::setEdgeRight(CircularEdgeSPtr edge_right) {
    this->edge_right_ = edge_right;
}

Plane3SPtr CircularArc::getSupportingPlane() const {
    DEBUG_SPTR(supporting_plane_);
    return supporting_plane_;
}

void CircularArc::setSupportingPlane(Plane3SPtr supporting_plane) {
    this->supporting_plane_ = supporting_plane;
}

SphericalSkeletonSPtr CircularArc::getSkel() const {
    DEBUG_WPTR(skel_);
    if (this->skel_.expired())
        return SphericalSkeletonSPtr();
    else
        return SphericalSkeletonSPtr(this->skel_);
}

void CircularArc::setSkel(SphericalSkeletonSPtr skel) {
    this->skel_ = skel;
}

std::list<CircularArcSPtr>::iterator CircularArc::getListIt() const {
    return this->list_it_;
}

void CircularArc::setListIt(std::list<CircularArcSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

bool CircularArc::hasNodeDst() const {
    bool result = false;
    if (node_dst_) {
        result = true;
    }
    return result;
}

std::string CircularArc::toString() const {
    std::string result("CircularArc(");
    result += util::StringFactory::fromPointer(this) + ", ";
    result += "src=" + node_src_->toString() + ", ";
    if (node_dst_) {
        result += "dst=" + node_dst_->toString();
    } else {
        result += "dir=<" + util::StringFactory::fromDouble((*direction_)[0]) + ", " +
                util::StringFactory::fromDouble((*direction_)[1]) + ", " +
                util::StringFactory::fromDouble((*direction_)[2]) + ">";
    }
    result += ")";
    return result;
}

} } }
