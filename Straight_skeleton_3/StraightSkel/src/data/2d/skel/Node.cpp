// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   data/2d/skel/Node.cpp
 * @author Gernot Walzl
 * @date   2012-02-03
 */

#include "data/2d/skel/Node.h"

#include "debug.h"
#include "data/2d/skel/Arc.h"
#include "util/StringFactory.h"

namespace data { namespace _2d { namespace skel {

Node::Node(Point2SPtr point) {
    height_ = 0.0;
    point_ = point;
    id_ = -1;
}

Node::~Node() {
    point_.reset();
}

NodeSPtr Node::create(Point2SPtr point) {
    NodeSPtr result = NodeSPtr(new Node(point));
    return result;
}

Point2SPtr Node::getPoint() const {
    CGAL_SS3_DEBUG_SPTR(point_);
    return point_;
}

void Node::setPoint(Point2SPtr point) {
    this->point_ = point;
}

CGAL::FT Node::getHeight() const {
    return height_;
}

void Node::setHeight(CGAL::FT height) {
    this->height_ = height;
}

StraightSkeletonSPtr Node::getSkel() const {
    return this->skel_.lock();
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

std::list<ArcWPtr>& Node::arcs() {
    return this->arcs_;
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

#ifdef USE_CGAL
CGAL::FT Node::getX() const { return this->point_->x(); }
CGAL::FT Node::getY() const { return this->point_->y(); }
#else
double Node::getX() const { return this->point_->getX(); }
double Node::getY() const { return this->point_->getY(); }
#endif

std::string Node::toString() const {
    std::string result("Node(");
    if (id_ != -1) {
        result += "id=" + util::StringFactory::fromInteger(id_) + ", ";
    } else {
        result += util::StringFactory::fromPointer(this) + ", ";
    }
    result += "<" + util::StringFactory::fromDouble(CGAL::to_double(getX())) + ", ";
    result += util::StringFactory::fromDouble(CGAL::to_double(getY())) + ">)";
    return result;
}

} } }
