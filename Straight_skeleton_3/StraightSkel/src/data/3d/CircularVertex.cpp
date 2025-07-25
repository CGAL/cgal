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
 * @file   data/3d/CircularVertex.cpp
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#include "data/3d/CircularVertex.h"

#include "debug.h"
#include "data/3d/CircularEdge.h"
#include "util/StringFactory.h"

namespace data { namespace _3d {

CircularVertex::CircularVertex(Point3SPtr point) {
    this->point_ = point;
    this->point_valid_ = true;
}

CircularVertex::~CircularVertex() {
    point_.reset();
}

CircularVertex::CircularVertex(const CircularVertex& vertex) {
    point_ = vertex.point_;
    point_valid_ = vertex.point_valid_;
}

CircularVertexSPtr CircularVertex::create(Point3SPtr point) {
    CircularVertexSPtr result = CircularVertexSPtr(new CircularVertex(point));
    return result;
}

CircularVertexSPtr CircularVertex::clone() const {
    CircularVertexSPtr result = CircularVertexSPtr(new CircularVertex(*this));
    return result;
}

Point3SPtr CircularVertex::getPoint() const {
    CGAL_SS3_DEBUG_SPTR(point_);
    return this->point_;
}

void CircularVertex::setPoint(Point3SPtr point) {
    this->point_ = point;
}

bool CircularVertex::isPointValid() const {
    return this->point_valid_;
}

void CircularVertex::setPointValid(bool point_valid) {
    this->point_valid_ = point_valid;
}

CircularEdgeSPtr CircularVertex::getEdgeIn() const {
    CGAL_SS3_DEBUG_WPTR(edge_in_);
    if (this->edge_in_.expired())
        return CircularEdgeSPtr();
    else
        return CircularEdgeSPtr(this->edge_in_);
}

void CircularVertex::setEdgeIn(CircularEdgeSPtr edge) {
    this->edge_in_ = edge;
}

CircularEdgeSPtr CircularVertex::getEdgeOut() const {
    CGAL_SS3_DEBUG_WPTR(edge_out_);
    if (this->edge_out_.expired())
        return CircularEdgeSPtr();
    else
        return CircularEdgeSPtr(this->edge_out_);
}

void CircularVertex::setEdgeOut(CircularEdgeSPtr edge) {
    this->edge_out_ = edge;
}

SphericalPolygonSPtr CircularVertex::getPolygon() const {
    CGAL_SS3_DEBUG_WPTR(polygon_);
    if (this->polygon_.expired())
        return SphericalPolygonSPtr();
    else
        return SphericalPolygonSPtr(this->polygon_);
}

void CircularVertex::setPolygon(SphericalPolygonSPtr polygon) {
    this->polygon_ = polygon;
}

std::list<CircularVertexSPtr>::iterator CircularVertex::getListIt() const {
    return this->list_it_;
}

void CircularVertex::setListIt(std::list<CircularVertexSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

CircularVertexDataSPtr CircularVertex::getData() const {
    CGAL_SS3_DEBUG_SPTR(data_);
    return this->data_;
}

void CircularVertex::setData(CircularVertexDataSPtr data) {
    this->data_ = data;
}

bool CircularVertex::hasData() const {
    bool result = false;
    if (data_) {
        result = true;
    }
    return result;
}

CircularVertexSPtr CircularVertex::next() const {
    CircularVertexSPtr result;
    if (!edge_out_.expired()) {
        CircularEdgeSPtr edge_out(edge_out_);
        result = edge_out->getVertexDst();
    }
    return result;
}

CircularVertexSPtr CircularVertex::prev() const {
    CircularVertexSPtr result;
    if (!edge_in_.expired()) {
        CircularEdgeSPtr edge_in(edge_in_);
        result = edge_in->getVertexSrc();
    }
    return result;
}

#ifdef USE_CGAL
CGAL::FT CircularVertex::getX() const { return this->point_->x(); }
CGAL::FT CircularVertex::getY() const { return this->point_->y(); }
CGAL::FT CircularVertex::getZ() const { return this->point_->z(); }
#else
double CircularVertex::getX() const { return this->point_->getX(); }
double CircularVertex::getY() const { return this->point_->getY(); }
double CircularVertex::getZ() const { return this->point_->getZ(); }
#endif

std::string CircularVertex::toString() const {
    std::string result("CircularVertex(");
    result += util::StringFactory::fromPointer(this) + ", ";
    result += "<" + util::StringFactory::fromDouble(CGAL::to_double(getX())) + ", ";
    result += util::StringFactory::fromDouble(CGAL::to_double(getY())) + ", ";
    result += util::StringFactory::fromDouble(CGAL::to_double(getZ())) + ">)";
    return result;
}

} }
