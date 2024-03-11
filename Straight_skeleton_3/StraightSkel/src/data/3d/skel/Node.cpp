/**
 * @file   data/3d/skel/Node.cpp
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#include "data/3d/skel/Node.h"

#include "data/3d/KernelFactory.h"
#include "data/3d/skel/Arc.h"
#include "debug.h"
#include "util/StringFactory.h"
#include <algorithm>

namespace data { namespace _3d { namespace skel {

Node::Node(Point3SPtr point) {
    offset_ = 0.0;
    point_ = point;
    id_ = -1;
}

Node::~Node() {
    point_.reset();
}

NodeSPtr Node::create(Point3SPtr point) {
    NodeSPtr result = NodeSPtr(new Node(point));
    return result;
}

Point3SPtr Node::getPoint() const {
    DEBUG_SPTR(point_);
    return point_;
}

void Node::setPoint(Point3SPtr point) {
    this->point_ = point;
}

double Node::getOffset() const {
    return offset_;
}

void Node::setOffset(double offset) {
    this->offset_ = offset;
}

StraightSkeletonSPtr Node::getSkel() const {
    DEBUG_WPTR(skel_);
    if (this->skel_.expired())
        return StraightSkeletonSPtr();
    else
        return StraightSkeletonSPtr(this->skel_);
}

void Node::setSkel(StraightSkeletonSPtr skel) {
    this->skel_ = skel;
}

std::list<NodeSPtr>::iterator Node::getListIt() const {
    return this->list_it_;
}

void Node::setListIt(std::list<NodeSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

int Node::getID() const {
    return this->id_;
}

void Node::setID(int id) {
    this->id_ = id;
}

void Node::addArc(ArcSPtr arc) {
    std::list<ArcWPtr>::iterator it = arcs_.insert(arcs_.end(), ArcWPtr(arc));
    if (arc->getNodeSrc() == shared_from_this()) {
        arc->setNodeSrcListIt(it);
    } else if (arc->hasNodeDst()) {
        if (arc->getNodeDst() == shared_from_this()) {
            arc->setNodeDstListIt(it);
        }
    }
}

bool Node::removeArc(ArcSPtr arc) {
    bool result = false;
    if (arc->getNodeSrc() == shared_from_this()) {
        arcs_.erase(arc->getNodeSrcListIt());
        arc->setNodeSrcListIt(std::list<ArcWPtr>::iterator());
        result = true;
    } else if (arc->hasNodeDst()) {
        if (arc->getNodeDst() == shared_from_this()) {
            arcs_.erase(arc->getNodeDstListIt());
            arc->setNodeDstListIt(std::list<ArcWPtr>::iterator());
            result = true;
        }
    }
    return result;
}

void Node::addSheet(SheetSPtr sheet) {
    sheets_.insert(sheets_.end(), SheetWPtr(sheet));
}

bool Node::removeSheet(SheetSPtr sheet) {
    bool result = false;
    std::list<SheetWPtr>::iterator it = sheets_.begin();
    while (it != sheets_.end()) {
        std::list<SheetWPtr>::iterator it_current = it;
        SheetWPtr sheet_wptr = *it++;
        if (!sheet_wptr.expired()) {
            if (sheet_wptr.lock() == sheet) {
                sheets_.erase(it_current);
                result = true;
                break;
            }
        }
    }
    return result;
}

bool Node::containsArc(ArcSPtr arc) const {
    ArcWPtr arc_wptr = ArcWPtr(arc);
    bool result = (arcs_.end() !=
            std::find(arcs_.begin(), arcs_.end(), arc_wptr));
    return result;
}

bool Node::containsSheet(SheetSPtr sheet) const {
    SheetWPtr sheet_wptr = SheetWPtr(sheet);
    bool result = (sheets_.end() !=
            std::find(sheets_.begin(), sheets_.end(), sheet_wptr));
    return result;
}

void Node::clear() {
    std::list<SheetWPtr>::iterator it_s = sheets_.begin();
    while (it_s != sheets_.end()) {
        SheetWPtr sheet_wptr = *it_s++;
        if (!sheet_wptr.expired()) {
            SheetSPtr sheet = SheetSPtr(sheet_wptr);
            removeSheet(sheet);
        }
    }
    sheets_.clear();
    std::list<ArcWPtr>::iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
        ArcWPtr arc_wptr = *it_a++;
        if (!arc_wptr.expired()) {
            ArcSPtr arc = ArcSPtr(arc_wptr);
            removeArc(arc);
        }
    }
    arcs_.clear();
}

std::list<ArcWPtr>& Node::arcs() {
    return this->arcs_;
}

std::list<SheetWPtr>& Node::sheets() {
    return this->sheets_;
}

unsigned int Node::degree() const {
    unsigned int result = 0;
    std::list<ArcWPtr>::const_iterator it_a = arcs_.begin();
    while (it_a != arcs_.end()) {
        ArcWPtr arc_wptr = *it_a++;
        if (!arc_wptr.expired()) {
            result++;
        }
    }
    return result;
}

double Node::getX() const {
#ifdef USE_CGAL
    return this->point_->x();
#else
    return this->point_->getX();
#endif
}

double Node::getY() const {
#ifdef USE_CGAL
    return this->point_->y();
#else
    return this->point_->getY();
#endif
}

double Node::getZ() const {
#ifdef USE_CGAL
    return this->point_->z();
#else
    return this->point_->getZ();
#endif
}

std::string Node::toString() const {
    std::string result("Node(");
    if (id_ != -1) {
        result += "id=" + util::StringFactory::fromInteger(id_) + ", ";
    } else {
        result += util::StringFactory::fromPointer(this) + ", ";
    }
    result += "<" + util::StringFactory::fromDouble(getX()) + ", ";
    result += util::StringFactory::fromDouble(getY()) + ", ";
    result += util::StringFactory::fromDouble(getZ()) + ">)";
    return result;
}

} } }
