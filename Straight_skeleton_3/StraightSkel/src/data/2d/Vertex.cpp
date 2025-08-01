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
 * @file   data/2d/Vertex.cpp
 * @author Gernot Walzl
 * @date   2011-11-22
 */

#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "debug.h"
#include "util/StringFactory.h"

#include <cmath>

namespace data { namespace _2d {

Vertex::Vertex(Point2SPtr point) {
    this->point_ = point;
    this->id_ = -1;
}

Vertex::~Vertex() {
    // intentionally does nothing
}

Vertex::Vertex(const Vertex& vertex) {
    point_ = vertex.point_;
    this->id_ = -1;
}

VertexSPtr Vertex::create(Point2SPtr point) {
    VertexSPtr result = VertexSPtr(new Vertex(point));
    return result;
}

VertexSPtr Vertex::clone() const {
    VertexSPtr result = VertexSPtr(new Vertex(*this));
    return result;
}

Point2SPtr Vertex::getPoint() const {
    CGAL_SS3_DEBUG_SPTR(point_);
    return this->point_;
}

void Vertex::setPoint(Point2SPtr point) {
    this->point_ = point;
}

EdgeSPtr Vertex::getEdgeIn() const {
    // CGAL_SS3_DEBUG_WPTR(edge_in_);
    if (this->edge_in_.expired())
        return EdgeSPtr();
    else
        return EdgeSPtr(this->edge_in_);
}

void Vertex::setEdgeIn(EdgeSPtr edge) {
    this->edge_in_ = edge;
}

EdgeSPtr Vertex::getEdgeOut() const {
    // CGAL_SS3_DEBUG_WPTR(edge_out_);
    if (this->edge_out_.expired())
        return EdgeSPtr();
    else
        return EdgeSPtr(this->edge_out_);
}

void Vertex::setEdgeOut(EdgeSPtr edge) {
    this->edge_out_ = edge;
}

PolygonSPtr Vertex::getPolygon() const {
    CGAL_SS3_DEBUG_WPTR(polygon_);
    if (this->polygon_.expired())
        return PolygonSPtr();
    else
        return PolygonSPtr(this->polygon_);
}

void Vertex::setPolygon(PolygonSPtr polygon) {
    this->polygon_ = polygon;
}

std::list<VertexSPtr>::iterator Vertex::getListIt() const {
    return this->list_it_;
}

void Vertex::setListIt(std::list<VertexSPtr>::iterator list_it) {
    this->list_it_ = list_it;
}

VertexDataSPtr Vertex::getData() const {
    CGAL_SS3_DEBUG_SPTR(data_);
    return this->data_;
}

void Vertex::setData(VertexDataSPtr data) {
    this->data_ = data;
}

bool Vertex::hasData() const {
    bool result = false;
    if (data_) {
        result = true;
    }
    return result;
}

int Vertex::getID() const {
    return this->id_;
}

void Vertex::setID(int id) {
    this->id_ = id;
}

VertexSPtr Vertex::next() const {
    VertexSPtr result;
    if (!edge_out_.expired()) {
        EdgeSPtr edge_out(edge_out_);
        result = edge_out->getVertexDst();
    }
    return result;
}

VertexSPtr Vertex::prev() const {
    VertexSPtr result;
    if (!edge_in_.expired()) {
        EdgeSPtr edge_in(edge_in_);
        result = edge_in->getVertexSrc();
    }
    return result;
}

#ifdef USE_CGAL
CGAL::FT Vertex::getX() const { return this->point_->x(); }
CGAL::FT Vertex::getY() const { return this->point_->y(); }
#else
double Vertex::getX() const { return this->point_->getX(); }
double Vertex::getY() const { return this->point_->getY(); }
#endif

double Vertex::angle() const {
    double result = 0.0;
    EdgeSPtr edge_in = getEdgeIn();
    EdgeSPtr edge_out = getEdgeOut();
    if (edge_in && edge_out) {
        VertexSPtr src = edge_in->getVertexSrc();
        VertexSPtr dst = edge_out->getVertexDst();
        if (src && dst) {
#ifdef USE_CGAL
            double arc_in = atan2(CGAL::to_double(src->getY() - getY()), CGAL::to_double(src->getX() - getX()));
            double arc_out = atan2(CGAL::to_double(dst->getY() - getY()), CGAL::to_double(dst->getX() - getX()));
#else
            double arc_in = atan2(src->getY() - getY(), src->getX() - getX());
            double arc_out = atan2(dst->getY() - getY(), dst->getX() - getX());
#endif
            result = arc_in - arc_out;
            if (result < 0.0) {
                result += 2*CGAL_PI;
            }
        }
    }
    return result;
}

bool Vertex::isReflex() const {
    bool result = (this->angle() > CGAL_PI);
    return result;
}

std::string Vertex::toString() const {
    std::string result("Vertex(");
    if (id_ != -1) {
        result += "id=" + util::StringFactory::fromInteger(id_) + ", ";
    } else {
        result += util::StringFactory::fromPointer(this) + ", ";
    }
    result += "<" + util::StringFactory::fromDouble(CGAL::to_double(getX())) + ", ";
    result += util::StringFactory::fromDouble(CGAL::to_double(getY())) + ">)";
    return result;
}

} }
