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
 * @file   data/3d/Facet.cpp
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#include "data/3d/Facet.h"

#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "data/2d/Polygon.h"
#include "data/2d/KernelFactory.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Triangle.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/skel/SkelFacetData.h"
#include "util/StringFactory.h"
#include "util/Configuration.h"

#include <CGAL/simplest_rational_in_interval.h>

#include <cmath>
#include <random>
#include <sstream>

namespace data { namespace _3d {

Facet::Facet() {
    this->id_ = -1;
}

Facet::~Facet() {
    edges_.clear();
    vertices_.clear();
}

FacetSPtr Facet::create() {
    FacetSPtr result = FacetSPtr(new Facet());
    return result;
}

FacetSPtr Facet::create(const std::vector<VertexSPtr>& vertices) {
    FacetSPtr result = FacetSPtr(new Facet());

    const std::size_t nv = vertices.size();
    for (std::size_t i=0; i<nv; ++i) {
        result->addVertex(vertices[i]);
    }
    for (std::size_t i=0; i<nv; ++i) {
        EdgeSPtr edge = vertices[i]->findEdge(vertices[(i+1)%nv]);
        if (!edge) {
            edge = Edge::create(vertices[i], vertices[(i+1)%nv]);
        }
        result->addEdge(edge);
    }
    return result;
}

FacetSPtr Facet::create(const std::vector<EdgeSPtr>& edges) {
    FacetSPtr result = FacetSPtr(new Facet());
    for (std::size_t i=0; i<edges.size(); ++i) {
        result->addEdge(edges[i]);
    }
    return result;
}

FacetSPtr Facet::clone() const {
    FacetSPtr result = Facet::create();
    std::map<VertexSPtr, VertexSPtr> old_to_new;
    std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        VertexSPtr vertex_c = vertex->clone();
        old_to_new[vertex] = vertex_c;
        result->addVertex(vertex_c);
    }
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr src = old_to_new.at(edge->getVertexSrc());
        VertexSPtr dst = old_to_new.at(edge->getVertexDst());
        EdgeSPtr edge_c = Edge::create(src, dst);
        CGAL_assertion(edge_c->getVertexSrc() == src && edge_c->getVertexDst() == dst);
        if (edge->getFacetL() == shared_from_this()) {
            edge_c->setFacetL(result);
        }
        if (edge->getFacetR() == shared_from_this()) {
            edge_c->setFacetR(result);
        }
        result->addEdge(edge_c);
    }

    result->plane_ = this->plane_;
    result->base_plane_ = this->base_plane_;
    result->final_plane_ = this->final_plane_;

    result->cachedSpeed_ = this->cachedSpeed_;
    result->cachedPlane_ = this->cachedPlane_;

    return result;
}

void Facet::addVertex(VertexSPtr vertex) {
    vertices_.insert(vertices_.end(), vertex);
    vertex->addFacet(shared_from_this());
}

bool Facet::removeVertex(VertexSPtr vertex) {
    bool result = false;
    std::list<VertexSPtr>::iterator it_v =
        std::find(vertices_.begin(), vertices_.end(), vertex);
    if (it_v != vertices_.end()) {
        vertices_.erase(it_v);
        vertex->removeFacet(shared_from_this());
        result = true;
    }
    return result;
}

bool Facet::hasVertex(VertexSPtr vertex) {
  std::list<VertexSPtr>::iterator it_v =
      std::find(vertices_.begin(), vertices_.end(), vertex);
  return (it_v != vertices_.end());
}

bool Facet::isTriangle() const {
  return (vertices_.size() == 3);
}

void Facet::addEdge(EdgeSPtr edge) {
    std::list<EdgeSPtr>::iterator it = edges_.insert(edges_.end(), edge);
    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    if (facet_l == shared_from_this() && facet_r == shared_from_this()) {
        std::list<EdgeSPtr>::iterator it_e =
            std::find(edges_.begin(), edges_.end(), edge);
        if (it_e == edge->getFacetLListIt()) {
            edge->setFacetRListIt(it);
        } else {
            edge->setFacetLListIt(it);
        }
    } else if (facet_l == shared_from_this()) {
        edge->setFacetLListIt(it);
    } else if (facet_r == shared_from_this()) {
        edge->setFacetRListIt(it);
    } else if (!facet_l) {
        edge->setFacetL(shared_from_this());
        edge->setFacetLListIt(it);
    } else if (!facet_r) {
        edge->setFacetR(shared_from_this());
        edge->setFacetRListIt(it);
    } else {
        CGAL_SS3_HDS_TRACE(edge->toString());
        throw std::runtime_error("The given edge already has a left and a right facet.");
    }
    VertexSPtr vertex_src = edge->src(shared_from_this());
    if (!containsVertex(vertex_src)) {
        addVertex(vertex_src);
    }
    VertexSPtr vertex_dst = edge->dst(shared_from_this());
    if (!containsVertex(vertex_dst)) {
        addVertex(vertex_dst);
    }
}

bool Facet::removeEdge(EdgeSPtr edge) {
    bool result = false;
    if (edge->getFacetL() == shared_from_this()) {
        edges_.erase(edge->getFacetLListIt());
        edge->setFacetL(FacetSPtr());
        edge->setFacetLListIt(std::list<EdgeSPtr>::iterator());
        result = true;
    } else if (edge->getFacetR() == shared_from_this()) {
        edges_.erase(edge->getFacetRListIt());
        edge->setFacetR(FacetSPtr());
        edge->setFacetRListIt(std::list<EdgeSPtr>::iterator());
        result = true;
    }
    return result;
}

EdgeSPtr Facet::findEdge(FacetSPtr facet) const {
    EdgeSPtr result = EdgeSPtr();
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
            result = edge;
            break;
        }
    }
    return result;
}

std::list<EdgeSPtr> Facet::findEdges(FacetSPtr facet) const {
    std::list<EdgeSPtr> result;
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getFacetL() == facet || edge->getFacetR() == facet) {
            result.push_back(edge);
        }
    }
    return result;
}

bool Facet::containsVertex(VertexSPtr vertex) const {
    bool result = (vertices_.end() !=
        std::find(vertices_.begin(), vertices_.end(), vertex));
    return result;
}

bool Facet::containsEdge(EdgeSPtr edge) const {
    bool result = (edges_.end() !=
        std::find(edges_.begin(), edges_.end(), edge));
    return result;
}

void Facet::sortVertices() {
    FacetSPtr self(shared_from_this());
    std::list<VertexSPtr> tmp;
    std::list<VertexSPtr>::iterator it_v_tmp;
    while (vertices_.size() > 0) {
        VertexSPtr first = VertexSPtr();
        VertexSPtr vertex = vertices_.front();
        EdgeSPtr edge = EdgeSPtr();
        std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
        while (it_e != edges_.end()) {
            EdgeSPtr my_edge = *it_e++;
            if (vertex == my_edge->src(self)) {
                edge = my_edge;
                break;
            }
        }
        while (vertex != first) {
            if (!first) {
                first = vertex;
                it_v_tmp = tmp.insert(tmp.end(), vertex);
            } else {
                tmp.insert(tmp.end(), vertex);
            }
            vertex = edge->dst(self);
            edge = edge->next(self);
        }
        while (it_v_tmp != tmp.end()) {
            VertexSPtr vertex = *it_v_tmp++;
            std::list<VertexSPtr>::iterator it_v =
                std::find(vertices_.begin(), vertices_.end(), vertex);
            if (it_v != vertices_.end()) {
                vertices_.erase(it_v);
            }
        }
    }
    vertices_.clear();
    it_v_tmp = tmp.begin();
    while (it_v_tmp != tmp.end()) {
        VertexSPtr vertex = *it_v_tmp++;
        vertices_.insert(vertices_.end(), vertex);
    }
}

void Facet::sortEdges() {
    std::list<EdgeSPtr> tmp;
    std::list<EdgeSPtr>::iterator it_e_tmp;
    while (edges_.size() > 0) {
        EdgeSPtr first = EdgeSPtr();
        EdgeSPtr edge = edges_.front();
        while (edge != first) {
            if (!first) {
                first = edge;
                it_e_tmp = tmp.insert(tmp.end(), edge);
            } else {
                tmp.insert(tmp.end(), edge);
            }
            edge = edge->next(shared_from_this());
        }
        while (it_e_tmp != tmp.end()) {
            EdgeSPtr edge = *it_e_tmp++;
            removeEdge(edge);
        }
    }
    edges_.clear();
    it_e_tmp = tmp.begin();
    while (it_e_tmp != tmp.end()) {
        EdgeSPtr edge = *it_e_tmp++;
        addEdge(edge);
    }
}

PolyhedronSPtr Facet::getPolyhedron() const {
    return this->polyhedron_.lock();
}

void Facet::setPolyhedron(PolyhedronSPtr polyhedron) {
    this->polyhedron_ = polyhedron;
}

std::list<FacetSPtr>::iterator Facet::getPolyhedronListIt() const {
    return this->polyhedron_list_it_;
}

void Facet::setPolyhedronListIt(std::list<FacetSPtr>::iterator list_it) {
    this->polyhedron_list_it_ = list_it;
}

FacetDataSPtr Facet::getData() const {
    CGAL_SS3_DEBUG_SPTR(this->data_);
    return this->data_;
}

void Facet::setData(FacetDataSPtr data) {
    this->data_ = data;
}

bool Facet::hasData() const {
    bool result = false;
    if (data_) {
        result = true;
    }
    return result;
}

std::list<VertexSPtr>& Facet::vertices() {
    return this->vertices_;
}

std::list<EdgeSPtr>& Facet::edges() {
    return this->edges_;
}

FacetSPtr Facet::next(VertexSPtr vertex) const {
    FacetSPtr result = FacetSPtr();
    std::list<FacetWPtr>::const_iterator it_f = vertex->facets().begin();
    while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f;
        if (FacetSPtr facet = facet_wptr.lock()) {
            if (facet == shared_from_this()) {
                if (vertex->degree() == 1) {
                    result = facet;
                }
                break;
            }
        }
        it_f++;
    }
    if (it_f != vertex->facets().end()) {
        std::list<FacetWPtr>::const_iterator it_f_begin = it_f++;
        if (it_f == vertex->facets().end()) {
            it_f = vertex->facets().begin();
        }
        while (it_f != it_f_begin) {
            FacetWPtr facet_wptr = *it_f++;
            if (it_f == vertex->facets().end()) {
                it_f = vertex->facets().begin();
            }
            if (FacetSPtr facet = facet_wptr.lock()) {
                std::list<EdgeWPtr>::const_iterator it_e = vertex->edges().begin();
                while (it_e != vertex->edges().end()) {
                    EdgeWPtr edge_wptr = *it_e++;
                    if (EdgeSPtr edge = edge_wptr.lock()) {
                        FacetSPtr facet_l = edge->getFacetL();
                        FacetSPtr facet_r = edge->getFacetR();
                        if ((facet_l == shared_from_this() &&
                                    facet_r == facet &&
                                    edge->getVertexDst() == vertex) ||
                                (facet_r == shared_from_this() &&
                                    facet_l == facet &&
                                    edge->getVertexSrc() == vertex)) {
                            result = facet;
                            break;
                        }
                    }
                }
            }
        }
    }
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
}

FacetSPtr Facet::prev(VertexSPtr vertex) const {
    FacetSPtr result = FacetSPtr();
    std::list<FacetWPtr>::const_reverse_iterator it_f = vertex->facets().rbegin();
    while (it_f != vertex->facets().rend()) {
        FacetWPtr facet_wptr = *it_f;
        if (FacetSPtr facet = facet_wptr.lock()) {
            if (facet == shared_from_this()) {
                if (vertex->degree() == 1) {
                    result = facet;
                }
                break;
            }
        }
        it_f++;
    }
    if (it_f != vertex->facets().rend()) {
        std::list<FacetWPtr>::const_reverse_iterator it_f_begin = it_f++;
        if (it_f == vertex->facets().rend()) {
            it_f = vertex->facets().rbegin();
        }
        while (it_f != it_f_begin) {
            FacetWPtr facet_wptr = *it_f++;
            if (it_f == vertex->facets().rend()) {
                it_f = vertex->facets().rbegin();
            }
            if (FacetSPtr facet = facet_wptr.lock()) {
                std::list<EdgeWPtr>::const_iterator it_e = vertex->edges().begin();
                while (it_e != vertex->edges().end()) {
                    EdgeWPtr edge_wptr = *it_e++;
                    if (EdgeSPtr edge = edge_wptr.lock()) {
                        FacetSPtr facet_l = edge->getFacetL();
                        FacetSPtr facet_r = edge->getFacetR();
                        if ((facet_l == shared_from_this() &&
                                    facet_r == facet &&
                                    edge->getVertexSrc() == vertex) ||
                                (facet_r == shared_from_this() &&
                                    facet_l == facet &&
                                    edge->getVertexDst() == vertex)) {
                            result = facet;
                            break;
                        }
                    }
                }
            }
        }
    }
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
}

void Facet::merge(FacetSPtr facet) {
    std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
    while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        facet->removeVertex(vertex);
        if (!containsVertex(vertex)) {
            addVertex(vertex);
        }
    }
    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if ((facet_l == facet && facet_r == shared_from_this()) ||
                (facet_r == facet && facet_l == shared_from_this())) {
            facet->removeEdge(edge);
            removeEdge(edge);
            edge->getPolyhedron()->removeEdge(edge);
        } else {
            if (facet_l == facet) {
                facet->removeEdge(edge);
                edge->setFacetL(shared_from_this());
                addEdge(edge);
            }
            if (facet_r == facet) {
                facet->removeEdge(edge);
                edge->setFacetR(shared_from_this());
                addEdge(edge);
            }
        }
    }

    CGAL_assertion(getPolyhedron()->isConsistent());
}

int Facet::getID() const {
    return this->id_;
}

void Facet::setID(int id) {
    this->id_ = id;
}

Plane3SPtr Facet::getPlane() const {
    CGAL_precondition(bool(this->plane_));
    return this->plane_;
}

void Facet::setPlane(Plane3SPtr plane) {
    this->plane_ = plane;
}

Plane3SPtr Facet::getBasePlane() const {
    CGAL_precondition(bool(this->base_plane_));
    return this->base_plane_;
}

void Facet::setBasePlane(Plane3SPtr plane) {
    this->base_plane_ = plane;
}

bool Facet::hasFinalPlane() const {
    return bool(this->final_plane_);
}

Plane3SPtr Facet::getFinalPlane() const {
    CGAL_precondition(bool(this->final_plane_));
    return this->final_plane_;
}

void Facet::setFinalPlane(Plane3SPtr plane) {
    this->final_plane_ = plane;
}

bool Facet::initPlane() {
    bool result = false;

    CGAL_SS3_HDS_TRACE("initPlane() of " << getID());

    Point3SPtr point_prev;
    std::vector<Point3SPtr> points;
    std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        Point3SPtr point = vertex->getPoint();
        if (point_prev != point) {
            points.push_back(point);
        }
        point_prev = point;
    }

    CGAL_SS3_HDS_TRACE("computing normals from " << points.size() << " points");
    CGAL_SS3_HDS_TRACE_CODE(for (std::size_t i=0; i<points.size(); ++i))
    CGAL_SS3_HDS_TRACE("point " << i << ": " << *points[i]);

    if (points.size() >= 3)
    {
        Point3SPtr p0 = points[0];
        Vector3SPtr normal = KernelFactory::createVector3(0,0,0);
        std::size_t last_i = points.size() - 1;
        for(std::size_t i=1; i<last_i; ++i) {
            Point3SPtr p1 = points[i];
            Point3SPtr p2 = points[i+1];
            CGAL_assertion(p1 && p2);
            *normal += 0.5 * CGAL::cross_product(*p2 - *p1, *p0 - *p1);
        }

        if(*normal != CGAL::NULL_VECTOR) {
            plane_ = KernelFactory::createPlane3(p0, normal);
            result = true;
        }
    }

    CGAL_postcondition(result);

    return result;
}

Plane3SPtr Facet::plane() {
    if (!this->plane_) {
        this->initPlane();
    }
    CGAL_SS3_DEBUG_SPTR(this->plane_);
    return this->plane_;
}

void Facet::normalizePlaneCoefficients()
{
    CGAL_precondition(bool(this->plane_));

    const CGAL::FT& a = plane_->a();
    const CGAL::FT& b = plane_->b();
    const CGAL::FT& c = plane_->c();
    const CGAL::FT& d = plane_->d();
    // this should be the only place with unavoidable SQRTs
    const CGAL::FT n = CGAL::approximate_sqrt(CGAL::square(a) + CGAL::square(b) + CGAL::square(c));

    if(!is_zero(n)) {
        // @todo to_double() it here too?
        plane_ = KernelFactory::createPlane3(a/n, b/n, c/n, d/n);
    }
}

bool Facet::isNormalizedPlane()
{
    CGAL_precondition(bool(this->plane_));
    const CGAL::FT& a = this->plane_->a();
    const CGAL::FT& b = this->plane_->b();
    const CGAL::FT& c = this->plane_->c();
    return (a*a + b*b + c*c - 1) <= 1e-5;
}

void Facet::storePlaneCoefficients()
{
    CGAL_precondition(bool(this->plane_));

    // need a different shared ptr here because plane_ will change with the perturbation
    cachedPlane_ = KernelFactory::createPlane3(*(plane_));

    CGAL_SS3_TRANSF_TRACE("caching plane of Facet " << this->id_ << " : [" << *cachedPlane_ << "]");
}

void Facet::perturbPlaneCoefficientsNudge(const double range) {
    CGAL_precondition(isNormalizedPlane());

    // @todo storing coefficients is only useful if we plan on untilting at the end.
    // storePlaneCoefficients();

    CGAL_SS3_TRANSF_TRACE("Nudging (Nudge) Face " << this->getID());
    CGAL_SS3_TRANSF_TRACE("  From coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                  << plane_->c() << " " << plane_->d() << "]");

    auto nudge = [&](const CGAL::FT& v) {
      static std::random_device rd;
      unsigned int s = 0; // rd()
      // CGAL_SS3_TRANSF_TRACE("seed = " << s);
      static std::mt19937 gen(s);
      static std::uniform_real_distribution<> rdist(-range, range);

      // Since we are perturbing, we might as well collapse the DAG of 'v'.
      // the point is also that once 'nv' is a double, its interval will be a singleton,
      // and we will have access to static filters
      double step = rdist(gen);
      double nv = CGAL::to_double(v) + step;
      return nv;
    };

    double na = nudge(plane_->a());
    double nb = nudge(plane_->b());
    double nc = nudge(plane_->c());
    double nd = nudge(plane_->d()); // @todo do not nudge 'd'? (mind the 'to_double()')

    double n = CGAL::approximate_sqrt(CGAL::square(na) + CGAL::square(nb) + CGAL::square(nc));
    CGAL_assertion(n != 0); // should not happen since we have normalized and the shift is tiny

    // below doesn't seem to matter? Probably need specific static filters...
#if 0
    plane_ = KernelFactory::createPlane3(na/n, nb/n, nc/n, nd/n);
#else
    // cast to_double() *after* the normalization to have double coordinates in the planes
    // the downside is that we won't have a^2 + b^2 + c^2 == 1,
    // but then again, who does...
    const double a = CGAL::to_double(na/n);
    const double b = CGAL::to_double(nb/n);
    const double c = CGAL::to_double(nc/n);
    const double d = CGAL::to_double(nd/n);
    plane_ = KernelFactory::createPlane3(a, b, c, d);
#endif

    CGAL_SS3_TRANSF_TRACE("  To coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                << plane_->c() << " " << plane_->d() << "]");

    CGAL_postcondition(isNormalizedPlane());
}

// I can go lower
void Facet::perturbPlaneCoefficientsSteps(int steps) {
    CGAL_precondition(isNormalizedPlane());

    // @todo storing coefficients is only useful if we plan on untilting at the end.
    // storePlaneCoefficients();

    CGAL_SS3_TRANSF_TRACE("Nudging (Steps) Face " << this->getID());
    CGAL_SS3_TRANSF_TRACE("  From coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                  << plane_->c() << " " << plane_->d() << "]");

    auto nudge = [&](const CGAL::FT& v) {
        std::random_device rd;
        unsigned int s = 0; // rd()
        // CGAL_SS3_TRANSF_TRACE("seed = " << s);
        std::mt19937 gen(s);
        std::uniform_int_distribution<> step_dis(1, steps);
        std::uniform_int_distribution<> dir_dis(0, 1); // 0: negative, 1: positive
        int n = step_dis(gen);
        int dir = dir_dis(gen);
        double direction = dir ? INFINITY : -INFINITY;
        double nv = CGAL::to_double(v);
        for (int i = 0; i < n; ++i) {
            nv = std::nextafter(nv, direction);
        }
        return nv;
    };

    double na = nudge(plane_->a());
    double nb = nudge(plane_->b());
    double nc = nudge(plane_->c());
    double nd = nudge(plane_->d()); // @todo do not nudge 'd'? (mind the 'to_double()')

    double n = CGAL::approximate_sqrt(CGAL::square(na) + CGAL::square(nb) + CGAL::square(nc));
    CGAL_assertion(n != 0); // should not happen since we have normalized and the shift is tiny

    // below doesn't seem to matter? Probably need specific static filters...
#if 0
    plane_ = KernelFactory::createPlane3(na/n, nb/n, nc/n, nd/n);
#else
    // cast to_double() *after* the normalization to have double coordinates in the planes
    // the downside is that we won't have a^2 + b^2 + c^2 == 1,
    // but then again, who does...
    plane_ = KernelFactory::createPlane3(CGAL::to_double(na/n),
                                         CGAL::to_double(nb/n),
                                         CGAL::to_double(nc/n),
                                         CGAL::to_double(nd/n));
#endif

    CGAL_SS3_TRANSF_TRACE("  To coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                << plane_->c() << " " << plane_->d() << "]");

    CGAL_postcondition(isNormalizedPlane());
}

// I can go lower
void Facet::perturbPlaneCoefficientsExact(const CGAL::FT& den) {
    CGAL_precondition(isNormalizedPlane());

    // @todo storing coefficients is only useful if we plan on untilting at the end.
    // storePlaneCoefficients();

    CGAL_SS3_TRANSF_TRACE("Nudging (Exact) Face " << this->getID());
    CGAL_SS3_TRANSF_TRACE("  From coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                  << plane_->c() << " " << plane_->d() << "]");

    auto nudge = [&](const CGAL::FT& v) {
        static std::random_device rd;
        unsigned int s = 0; // rd()
        // CGAL_SS3_TRANSF_TRACE("seed = " << s);
        static std::mt19937 gen(s);
        static std::uniform_real_distribution<> rdist(0, 1);
        CGAL::FT step = rdist(gen);
        return v + step / CGAL::square(den);
    };

    CGAL::FT na = nudge(plane_->a());
    CGAL::FT nb = nudge(plane_->b());
    CGAL::FT nc = nudge(plane_->c());
    CGAL::FT nd = nudge(plane_->d()); // @todo do not nudge 'd'?

    // so small, needless to normalize (@todo should we even normalize others?)
    plane_ = KernelFactory::createPlane3(na, nb, nc, nd);

    CGAL_SS3_TRANSF_TRACE("  To coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                << plane_->c() << " " << plane_->d() << "]");

    CGAL_postcondition(isNormalizedPlane());
}

void Facet::perturbPlaneCoefficientsFixedPoints(const double range,
                                                const std::vector<Point3SPtr>& fixed_points) {
    CGAL_precondition(isNormalizedPlane());
    CGAL_precondition(fixed_points.size() <= 2);

    // @todo storing coefficients is only useful if we plan on untilting at the end.
    // storePlaneCoefficients();

    CGAL_SS3_TRANSF_TRACE("Nudging Face " << this->getID());
    CGAL_SS3_TRANSF_TRACE("  From coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                  << plane_->c() << " " << plane_->d() << "]");
    CGAL_SS3_TRANSF_TRACE("  with " << fixed_points.size() << " fixed points");
    CGAL_SS3_TRANSF_TRACE_CODE(for (Point3SPtr fp : fixed_points))
    CGAL_SS3_TRANSF_TRACE("    " << *fp);

    static std::random_device rd;
    unsigned int s = 0; // rd()
    // CGAL_SS3_TRANSF_TRACE("seed = " << s);
    static std::mt19937 gen(s);
    static std::uniform_real_distribution<> rdist(-range, range);

    auto nudge = [&](const CGAL::FT& v) {
        // Since we are perturbing, we might as well collapse the DAG of 'v'.
        // the point is also that once 'nv' is a double, its interval will be a singleton,
        // and we will have access to static filters
        double step = rdist(gen);
        double nv = CGAL::to_double(v) + step;
        return nv;
    };

#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
    auto nudge_to_simplest_rational_in_interval = [&](const CGAL::FT& v) {
        double d1 = nudge(v);
        double d2 = nudge(v);
        if (d2 < d1) {
            std::swap(d1, d2);
        }
        CGAL::FT nv = CGAL::simplest_rational_in_interval<CGAL::K::Exact_kernel::FT>(d1, d2);
        return nv;
    };
#endif

    // 0 fixed points: nudge all coefficients independently
    if (fixed_points.size() == 0) {
#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        CGAL::FT na = nudge_to_simplest_rational_in_interval(plane_->a());
        CGAL::FT nb = nudge_to_simplest_rational_in_interval(plane_->b());
        CGAL::FT nc = nudge_to_simplest_rational_in_interval(plane_->c());
        CGAL::FT nd = nudge_to_simplest_rational_in_interval(plane_->d());
#else
        double na = nudge(plane_->a());
        double nb = nudge(plane_->b());
        double nc = nudge(plane_->c());
        double nd = nudge(plane_->d());
#endif
        plane_ = KernelFactory::createPlane3(na, nb, nc, nd);
    } else if (fixed_points.size() == 1) {
        // 1 fixed point: nudge (a, b, c), recompute d so the plane passes through the point
        Point3SPtr p0 = fixed_points[0];

#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        CGAL::FT na = nudge_to_simplest_rational_in_interval(plane_->a());
        CGAL::FT nb = nudge_to_simplest_rational_in_interval(plane_->b());
        CGAL::FT nc = nudge_to_simplest_rational_in_interval(plane_->c());
#else
        double na = nudge(plane_->a());
        double nb = nudge(plane_->b());
        double nc = nudge(plane_->c());
#endif
        const CGAL::FT& x0 = p0->x();
        const CGAL::FT& y0 = p0->y();
        const CGAL::FT& z0 = p0->z();
        CGAL::FT d = - (na * x0 + nb * y0 + nc * z0);
        plane_ = KernelFactory::createPlane3(na, nb, nc, d);
        CGAL_postcondition(plane_->has_on(*p0));
    } else if (fixed_points.size() == 2) {
        // 2 fixed points: construct a plane through both points, nudge the normal within the allowed family
        Point3SPtr p0 = fixed_points[0];
        Point3SPtr p1 = fixed_points[1];
        CGAL_assertion(*p0 != *p1);

        const CGAL::FT& p0x = p0->x();
        const CGAL::FT& p0y = p0->y();
        const CGAL::FT& p0z = p0->z();
        const CGAL::FT& p1x = p1->x();
        const CGAL::FT& p1y = p1->y();
        const CGAL::FT& p1z = p1->z();

        // Step 1: Direction vector between points
        CGAL::FT ux = p1x - p0x;
        CGAL::FT uy = p1y - p0y;
        CGAL::FT uz = p1z - p0z;
        CGAL::FT uu = ux*ux + uy*uy + uz*uz;

        // Step 2: Original normal
        const CGAL::FT& a0 = plane_->a();
        const CGAL::FT& b0 = plane_->b();
        const CGAL::FT& c0 = plane_->c();

        // Step 3: Project original normal onto plane orthogonal to u
        CGAL::FT dot = a0*ux + b0*uy + c0*uz;
        CGAL::FT ab = a0 - dot * ux / uu;
        CGAL::FT bb = b0 - dot * uy / uu;
        CGAL::FT cb = c0 - dot * uz / uu;

        // Step 4: Find a direction to nudge (cross product)
        CGAL::FT vx = uy * cb - uz * bb;
        CGAL::FT vy = uz * ab - ux * cb;
        CGAL::FT vz = ux * bb - uy * ab;

        // Step 5: Nudge the normal
#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        CGAL::FT epsilon = nudge_to_simplest_rational_in_interval(rdist(gen));
#else
        double epsilon = rdist(gen);
#endif

        CGAL::FT a1 = ab + epsilon * vx;
        CGAL::FT b1 = bb + epsilon * vy;
        CGAL::FT c1 = cb + epsilon * vz;

        // Step 6: Compute d so plane passes through p0
        CGAL::FT d1 = - (a1 * p0x + b1 * p0y + c1 * p0z);
        plane_ = KernelFactory::createPlane3(a1, b1, c1, d1);

        CGAL_postcondition(plane_->has_on(*p0));
        CGAL_postcondition(plane_->has_on(*p1));
    } else {
        CGAL_SS3_TRANSF_TRACE("Error: called fixed point facet perturbation with > 2 fixed points");
    }

    CGAL_SS3_TRANSF_TRACE("  To coefficients [" << plane_->a() << " " << plane_->b() << " "
                                                << plane_->c() << " " << plane_->d() << "]");

    CGAL_postcondition(isNormalizedPlane());
}

void Facet::perturbPlaneCoefficientsHighDegrees(const double range) {
    std::vector<Point3SPtr> high_degree_points;
    for (VertexSPtr v : vertices_) {
        if (v->degree() > 3) {
            high_degree_points.push_back(v->getPoint());
        }
    }
    return perturbPlaneCoefficientsFixedPoints(range, high_degree_points);
}

void Facet::perturbPlaneCoefficients(PerturbationType type)
{
    double range = 1e-10;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        range = config->getDouble("main", "rand_move_points_range");
    }

    if (type == PerturbationType::NUDGE) {
        perturbPlaneCoefficientsNudge(range);
    } else if (type == PerturbationType::STEPS) {
        perturbPlaneCoefficientsSteps(10);
    } else if (type == PerturbationType::EXACT) {
        perturbPlaneCoefficientsExact(CGAL::FT(1e100));
    } else if (type == PerturbationType::HIGH_DEGREES) {
        perturbPlaneCoefficientsHighDegrees(range);
    } else {
        CGAL_error_msg("Unknown perturbation type");
    }
}

void Facet::restorePlaneCoefficients(const CGAL::FT& perturbationOffset,
                                     const CGAL::FT& perturbationEndOffset)
{
    CGAL_SS3_TRANSF_TRACE("plane of Facet " << this->id_ << " is [" << *plane_ << "]");

    if (!cachedPlane_) {
        std::cerr << "Warning: no plane coefficients to restore" << std::endl;
        return;
    }

    CGAL::FT speed = 1.0;
    if (hasData()) {
        speed = std::dynamic_pointer_cast<data::_3d::skel::SkelFacetData>(getData())->getSpeed();
    }

    CGAL_SS3_TRANSF_TRACE("OLD d = " << cachedPlane_->d());
    CGAL_SS3_TRANSF_TRACE("perturbationOffset = " << perturbationOffset);
    CGAL_SS3_TRANSF_TRACE("perturbationEndOffset = " << perturbationEndOffset);

    // The minus sign "d - ..." is because we shrink, so the plane needs to be offset
    // by the difference of offsets, but in the direction opposite of its normal.
    //
    // This is similar to when we call, e.g.:
    //   Plane3SPtr offset_plane_l = KernelWrapper::offsetPlane(plane_l, - speed_l);
    //                                                                  ^^^
    CGAL:: FT d = cachedPlane_->d() - speed * (perturbationEndOffset - perturbationOffset);

    plane_ = KernelFactory::createPlane3(cachedPlane_->a(), cachedPlane_->b(), cachedPlane_->c(), d);
    CGAL_assertion_code(CGAL::FT sq_n = CGAL::square(plane_->a()) + CGAL::square(plane_->b()) + CGAL::square(plane_->c()));
    CGAL_assertion((sq_n - 1) < 1e-5);

    CGAL_SS3_TRANSF_TRACE("plane of Facet " << this->id_ << " restored to [" << *plane_ << "]");
}

bool Facet::makeFirstConvex() {
    bool result = false;
    if (!plane_) {
        return false;
    }
    EdgeSPtr edge_begin;
    FacetSPtr self(shared_from_this());
    Vector3SPtr normal = KernelFactory::createVector3(plane_);
    EdgeSPtr edge = edges_.front();
    EdgeSPtr first = EdgeSPtr();
    while (edge != first) {
        if (!first) {
            first = edge;
        }
        EdgeSPtr edge_next = edge->next(self);
        Point3SPtr points[3];
        points[0] = edge->src(self)->getPoint();
        points[1] = edge->dst(self)->getPoint();
        points[2] = edge_next->dst(self)->getPoint();

#if 1 // @fixme is this correct? the CGAL_PI/4.0 below is confusing... Was it supposed to be CGAL_PI/2.0?
        if(!CGAL::collinear(*(points[0]), *(points[1]), *(points[2])))
        {
          if(CGAL::angle(*(points[0]), *(points[1]), *(points[2]), *normal) == CGAL::ACUTE)
          {
            edge_begin = edge;
            result = true;
            break;
          }
        }
#else // old code; has issues with collinear points
        if (points[0] != points[1] &&
                points[1] != points[2] &&
                points[2] != points[0])
        {
          Plane3SPtr plane_current = KernelFactory::createPlane3(
                  points[0], points[1], points[2]);
          Vector3SPtr normal_current = KernelFactory::createVector3(plane_current);
          double angle = 0.0;
          double arg = 0.0;
          arg = CGAL::to_double(((*normal)*(*normal_current)) /
                  CGAL::disallowed_sqrt(normal->squared_length() * normal_current->squared_length()));

          // fixes issues with floating point precision
          if (arg <= -1.0) {
              angle = CGAL_PI;
          } else if (arg >= 1.0) {
              angle = 0.0;
          } else {
              angle = acos(arg);
          }
          if (angle < CGAL_PI/4.0) {
              edge_begin = edge;
              result = true;
              break;
          }
        }
#endif
        edge = edge_next;
    }

    if (edge_begin) {
        std::list<EdgeSPtr>::iterator it_e = edges_.insert(edges_.begin(), edge_begin);
        if (edge_begin->getFacetL() == shared_from_this()) {
            edges_.erase(edge->getFacetLListIt());
            edge_begin->setFacetLListIt(it_e);
        } else if (edge_begin->getFacetR() == shared_from_this()) {
            edges_.erase(edge->getFacetRListIt());
            edge_begin->setFacetRListIt(it_e);
        }
        sortEdges();
        VertexSPtr vertex_begin = edge_begin->src(shared_from_this());
        std::list<VertexSPtr>::iterator it_v =
            std::find(vertices_.begin(), vertices_.end(), vertex_begin);
        vertices_.erase(it_v);
        vertices_.insert(vertices_.begin(), vertex_begin);
        sortVertices();
    }
    if (!result) {
        CGAL_SS3_HDS_TRACE("Warning: Unable to make first 3 vertices convex.");
        CGAL_SS3_HDS_TRACE(toString());
    }

    return result;
}

int Facet::numHighDegreeVertices() const {
    int result = 0;
    // get min degree
    for (VertexSPtr v : vertices_) {
        if (v->degree() > 3) {
            ++result;
        }
    }
    return result;
}

/**
 * http://de.wikipedia.org/wiki/Drehmatrix
 *
 *              [ n_x^2 (1 - cos(alpha)) + cos(alpha)         n_x n_y (1 - cos(alpha)) - n_z sin(alpha)   n_x n_z (1 - cos(alpha)) + n_y sin(alpha) ]
 * R_n(alpha) = [ n_y n_x (1 - cos(alpha)) + n_z sin(alpha)   n_y^2 (1 - cos(alpha)) + cos(alpha)         n_y n_z (1 - cos(alpha)) - n_x sin(alpha) ]
 *              [ n_z n_x (1 - cos(alpha)) - n_y sin(alpha)   n_z n_y (1 - cos(alpha)) + n_x sin(alpha)   n_z^2 (1 - cos(alpha)) + cos(alpha)       ]
 *
 * n_z = 0
 */
data::_2d::PolygonSPtr Facet::toPolygon() {
    data::_2d::PolygonSPtr result = data::_2d::Polygon::create();
    Vector3SPtr normal = KernelFactory::createVector3(plane());
    CGAL::FT n[3];   // n = normal.cross([0 0 1])
    n[0] = (*normal)[1];
    n[1] = -(*normal)[0];
    n[2] = 0.0;
    CGAL::FT length = 0.0;
    if ((*normal)[0] != 0.0 || (*normal)[1] != 0.0) {
        // normalize n
        for (unsigned int i = 0; i < 3; i++) {
            length += n[i]*n[i];
        }
        length = CGAL::disallowed_sqrt(length);
        for (unsigned int i = 0; i < 3; i++) {
            n[i] /= length;
        }
    }
    length = 0.0;
    for (unsigned int i = 0; i < 3; i++) {
        length += (*normal)[i] * (*normal)[i];
    }
    length = CGAL::disallowed_sqrt(length);
    // alpha = acos((normal * [0 0 1]) / normal.length())
    double alpha = acos(CGAL::to_double((*normal)[2] / length));
    const CGAL::FT cos_alpha = std::cos(alpha);
    const CGAL::FT sin_alpha = std::sin(alpha);
    Vector3SPtr r_x = KernelFactory::createVector3(
            n[0] * n[0] * (1 - cos_alpha) + cos_alpha,
            n[0] * n[1] * (1 - cos_alpha) - n[2] * sin_alpha,
            n[0] * n[2] * (1 - cos_alpha) + n[1] * sin_alpha);
    Vector3SPtr r_y = KernelFactory::createVector3(
            n[1] * n[0] * (1 - cos_alpha) + n[2] * sin_alpha,
            n[1] * n[1] * (1 - cos_alpha) + cos_alpha,
            n[1] * n[2] * (1 - cos_alpha) - n[0] * sin_alpha);
    std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        Vector3SPtr v_p = KernelFactory::createVector3(vertex->getPoint());
        CGAL::FT x2 = (*r_x) * (*v_p);
        CGAL::FT y2 = (*r_y) * (*v_p);
        data::_2d::Point2SPtr p2 = data::_2d::KernelFactory::createPoint2(x2, y2);
        data::_2d::VertexSPtr vertex2 = data::_2d::Vertex::create(p2);
        result->addVertex(vertex2);
    }
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr src = edge->src(shared_from_this());
        VertexSPtr dst = edge->dst(shared_from_this());
        data::_2d::VertexSPtr src2;
        data::_2d::VertexSPtr dst2;
        std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
        std::list<data::_2d::VertexSPtr>::const_iterator it_v2 = result->vertices().begin();
        while (it_v != vertices_.end() && it_v2 != result->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            data::_2d::VertexSPtr vertex2 = *it_v2++;
            if (src == vertex) {
                src2 = vertex2;
            }
            if (dst == vertex) {
                dst2 = vertex2;
            }
        }
        data::_2d::EdgeSPtr edge2 = data::_2d::Edge::create(src2, dst2);
        result->addEdge(edge2);
    }
    return result;
}

std::string Facet::toString() const {
    std::stringstream sstr;
    sstr << "Facet(";
    sstr << "id=" << util::StringFactory::fromInteger(id_) << ", ";
    // sstr << "addr=" << util::StringFactory::fromPointer(this) << ", ";
    if (plane_) {
        sstr << "Plane: <" << util::StringFactory::fromDouble(CGAL::to_double(plane_->a())) << ", "
             << util::StringFactory::fromDouble(CGAL::to_double(plane_->b())) << ", "
             << util::StringFactory::fromDouble(CGAL::to_double(plane_->c())) << ", "
             << util::StringFactory::fromDouble(CGAL::to_double(plane_->d())) << ">, ";
    }

    if (hasData()) {
        sstr << "Speed: " << std::dynamic_pointer_cast<data::_3d::skel::SkelFacetData>(getData())->getSpeed() << ", ";
    }

    sstr << "Vertices:" + util::StringFactory::fromInteger(vertices_.size()) + ", ";
    sstr << "Edges:" + util::StringFactory::fromInteger(edges_.size()) + ",";
    if (vertices_.size() > 0) {
        sstr << "\n";
        std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
        while (it_v != vertices_.end()) {
            VertexSPtr vertex = *it_v++;
            sstr << vertex->toString() << "\n";
        }
    }
    if (edges_.size() > 0) {
        sstr << "\n";
        std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
            while (it_e != edges_.end()) {
            EdgeSPtr edge = *it_e++;
            sstr << edge->toString() << "\n";
        }
    }
    sstr << ") END FACET TOSTRING()\n";

    return sstr.str();
}

} }
