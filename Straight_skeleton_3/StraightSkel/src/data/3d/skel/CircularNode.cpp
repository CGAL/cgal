/**
 * @file   data/3d/skel/CircularNode.cpp
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#include "data/3d/skel/CircularNode.h"

#include "debug.h"
#include "data/3d/skel/CircularArc.h"
#include "util/StringFactory.h"
#include <algorithm>

namespace data { namespace _3d { namespace skel {

CircularNode::CircularNode(Point3SPtr point) {
    offset_ = 0.0;
    point_ = point;
}

CircularNode::~CircularNode() {
    point_.reset();
    skel_.reset();
}

CircularNodeSPtr CircularNode::create(Point3SPtr point) {
    CircularNodeSPtr result = CircularNodeSPtr(new CircularNode(point));
    return result;
}

Point3SPtr CircularNode::getPoint() const {
    DEBUG_SPTR(point_);
    return point_;
}

void CircularNode::setPoint(Point3SPtr point) {
    this->point_ = point;
}

double CircularNode::getOffset() const {
    return offset_;
}

void CircularNode::setOffset(double offset) {
    this->offset_ = offset;
}

SphericalSkeletonSPtr CircularNode::getSkel() const {
    DEBUG_WPTR(skel_);
    if (this->skel_.expired())
        return SphericalSkeletonSPtr();
    else
        return SphericalSkeletonSPtr(this->skel_);
}

void CircularNode::setSkel(SphericalSkeletonSPtr skel) {
    this->skel_ = skel;
}

std::list<CircularNodeSPtr>::iterator CircularNode::getListIt() const {
    return this->list_it_;
}

void CircularNode::setListIt(std::list<CircularNodeSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

void CircularNode::addArc(CircularArcSPtr arc) {
    std::list<CircularArcWPtr>::iterator it = arcs_.insert(arcs_.end(), CircularArcWPtr(arc));
    if (arc->getNodeSrc() == shared_from_this()) {
        arc->setNodeSrcListIt(it);
    } else if (arc->hasNodeDst()) {
        if (arc->getNodeDst() == shared_from_this()) {
            arc->setNodeDstListIt(it);
        }
    }
}

bool CircularNode::removeArc(CircularArcSPtr arc) {
    bool result = false;
    if (arc->getNodeSrc() == shared_from_this()) {
        arcs_.erase(arc->getNodeSrcListIt());
        arc->setNodeSrcListIt(std::list<CircularArcWPtr>::iterator());
        result = true;
    } else if (arc->hasNodeDst()) {
        if (arc->getNodeDst() == shared_from_this()) {
            arcs_.erase(arc->getNodeDstListIt());
            arc->setNodeDstListIt(std::list<CircularArcWPtr>::iterator());
            result = true;
        }
    }
    return result;
}

bool CircularNode::containsArc(CircularArcSPtr arc) const {
    CircularArcWPtr arc_wptr = CircularArcWPtr(arc);
    bool result = (arcs_.end() !=
            std::find(arcs_.begin(), arcs_.end(), arc_wptr));
    return result;
}

void CircularNode::clear() {
    std::list<CircularArcWPtr>::iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
        CircularArcWPtr arc_wptr = *it_a++;
        if (!arc_wptr.expired()) {
            CircularArcSPtr arc = CircularArcSPtr(arc_wptr);
            removeArc(arc);
        }
    }
    arcs_.clear();
}

std::list<CircularArcWPtr>& CircularNode::arcs() {
    return this->arcs_;
}

unsigned int CircularNode::degree() {
    unsigned int result = 0;
    std::list<CircularArcWPtr>::iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
        CircularArcWPtr arc_wptr = *it_a++;
        if (!arc_wptr.expired()) {
            result++;
        }
    }
    return result;
}

double CircularNode::getX() const {
#ifdef USE_CGAL
    return this->point_->x();
#else
    return this->point_->getX();
#endif
}

double CircularNode::getY() const {
#ifdef USE_CGAL
    return this->point_->y();
#else
    return this->point_->getY();
#endif
}

double CircularNode::getZ() const {
#ifdef USE_CGAL
    return this->point_->z();
#else
    return this->point_->getZ();
#endif
}

std::string CircularNode::toString() const {
    std::string result("CircularNode(");
    result += util::StringFactory::fromPointer(this) + ", ";
    result += "<" + util::StringFactory::fromDouble(getX()) + ", ";
    result += util::StringFactory::fromDouble(getY()) + ", ";
    result += util::StringFactory::fromDouble(getZ()) + ">)";
    return result;
}

} } }
