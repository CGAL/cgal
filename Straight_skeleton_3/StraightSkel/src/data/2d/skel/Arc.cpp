/**
 * @file   data/2d/skel/Arc.cpp
 * @author Gernot Walzl
 * @date   2012-02-02
 */

#include "data/2d/skel/Arc.h"

#include "debug.h"
#include "data/2d/KernelFactory.h"
#include "data/2d/skel/Node.h"
#include "util/StringFactory.h"

namespace data { namespace _2d { namespace skel {

Arc::Arc(NodeSPtr node_src, Vector2SPtr direction) {
    node_src_ = node_src;
    direction_ = direction;
    id_ = -1;
}

Arc::Arc(NodeSPtr node_src, NodeSPtr node_dst) {
    node_src_ = node_src;
    node_dst_ = node_dst;
    id_ = -1;
}

Arc::~Arc() {
    node_src_.reset();
    node_dst_.reset();
    direction_.reset();
    edge_left_.reset();
    edge_right_.reset();
}

ArcSPtr Arc::create(NodeSPtr node_src, Vector2SPtr direction) {
    ArcSPtr result = ArcSPtr(new Arc(node_src, direction));
    node_src->addArc(result);
    return result;
}

ArcSPtr Arc::create(NodeSPtr node_src, NodeSPtr node_dst) {
    ArcSPtr result = ArcSPtr(new Arc(node_src, node_dst));
    node_src->addArc(result);
    node_dst->addArc(result);
    return result;
}

NodeSPtr Arc::getNodeSrc() const {
    DEBUG_SPTR(node_src_);
    return node_src_;
}

void Arc::setNodeSrc(NodeSPtr node_src) {
    this->node_src_ = node_src;
}

std::list<ArcWPtr>::iterator Arc::getNodeSrcListIt() const {
    return this->node_src_list_it_;
}

void Arc::setNodeSrcListIt(std::list<ArcWPtr>::iterator node_src_list_it) {
    this->node_src_list_it_ = node_src_list_it;
}

NodeSPtr Arc::getNodeDst() const {
    DEBUG_SPTR(node_dst_);
    return node_dst_;
}

void Arc::setNodeDst(NodeSPtr node_dst) {
    this->node_dst_ = node_dst;
}

std::list<ArcWPtr>::iterator Arc::getNodeDstListIt() const {
    return this->node_dst_list_it_;
}

void Arc::setNodeDstListIt(std::list<ArcWPtr>::iterator node_dst_list_it) {
    this->node_dst_list_it_ = node_dst_list_it;
}

Vector2SPtr Arc::getDirection() const {
    DEBUG_SPTR(direction_);
    return this->direction_;
}

void Arc::setDirection(Vector2SPtr direction) {
    this->direction_ = direction;
}

EdgeSPtr Arc::getEdgeLeft() const {
    DEBUG_SPTR(edge_left_);
    return edge_left_;
}

void Arc::setEdgeLeft(EdgeSPtr edge_left) {
    this->edge_left_ = edge_left;
}

EdgeSPtr Arc::getEdgeRight() const {
    DEBUG_SPTR(edge_right_);
    return edge_right_;
}

void Arc::setEdgeRight(EdgeSPtr edge_right) {
    this->edge_right_ = edge_right;
}

StraightSkeletonSPtr Arc::getSkel() const {
    DEBUG_WPTR(skel_);
    if (this->skel_.expired())
        return StraightSkeletonSPtr();
    else
        return StraightSkeletonSPtr(this->skel_);
}

void Arc::setSkel(StraightSkeletonSPtr skel) {
    this->skel_ = skel;
}

std::list<ArcSPtr>::iterator Arc::getListIt() const {
    return this->list_it_;
}

void Arc::setListIt(std::list<ArcSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

int Arc::getID() const {
    return this->id_;
}

void Arc::setID(int id) {
    this->id_ = id;
}

Line2SPtr Arc::line() const {
    Line2SPtr result = Line2SPtr();
    if (node_dst_) {
        result = KernelFactory::createLine2(
            node_src_->getPoint(), node_dst_->getPoint());
    } else if (direction_) {
        result = KernelFactory::createLine2(
            node_src_->getPoint(), direction_);
    }
    return result;
}

bool Arc::hasNodeDst() const {
    bool result = false;
    if (node_dst_) {
        result = true;
    }
    return result;
}

std::string Arc::toString() const {
    std::string result("Arc(");
    if (id_ != -1) {
        result += "id=" + util::StringFactory::fromInteger(id_) + ", ";
    } else {
        result += util::StringFactory::fromPointer(this) + ", ";
    }
    result += "src=" + node_src_->toString() + ", ";
    if (node_dst_) {
        result += "dst=" + node_dst_->toString();
    } else {
        result += "dir=<" + util::StringFactory::fromDouble((*direction_)[0]) + ", " +
                util::StringFactory::fromDouble((*direction_)[1]) + ">";
    }
    result += ")";
    return result;
}

} } }
