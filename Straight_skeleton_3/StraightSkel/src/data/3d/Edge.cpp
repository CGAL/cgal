/**
 * @file   data/3d/Edge.cpp
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#include "data/3d/Edge.h"

#include "debug.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/Vertex.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "util/StringFactory.h"
#include <cmath>

namespace data { namespace _3d {

Edge::Edge(VertexSPtr src, VertexSPtr dst) {
    vertex_src_ = src;
    vertex_dst_ = dst;
    id_ = -1;
}

Edge::~Edge() {
    vertex_src_.reset();
    vertex_dst_.reset();
}

Edge::Edge(const Edge& edge) {
    vertex_src_ = edge.vertex_src_->clone();
    vertex_dst_ = edge.vertex_dst_->clone();
    id_ = -1;
}

EdgeSPtr Edge::create(VertexSPtr src, VertexSPtr dst) {
    EdgeSPtr result = EdgeSPtr(new Edge(src, dst));
    src->addEdge(result);
    dst->addEdge(result);
    return result;
}

EdgeSPtr Edge::clone() const {
    EdgeSPtr result = EdgeSPtr(new Edge(*this));
    result->vertex_src_->addEdge(result);
    result->vertex_dst_->addEdge(result);
    return result;
}

VertexSPtr Edge::getVertexSrc() const {
    DEBUG_SPTR(this->vertex_src_);
    return this->vertex_src_;
}

void Edge::setVertexSrc(VertexSPtr src) {
    this->vertex_src_ = src;
}

std::list<EdgeWPtr>::iterator Edge::getVertexSrcListIt() const {
    return this->vertex_src_list_it_;
}

void Edge::setVertexSrcListIt(std::list<EdgeWPtr>::iterator list_it) {
    this->vertex_src_list_it_ = list_it;
}

VertexSPtr Edge::getVertexDst() const {
    DEBUG_SPTR(this->vertex_dst_);
    return this->vertex_dst_;
}

void Edge::setVertexDst(VertexSPtr dst) {
    this->vertex_dst_ = dst;
}

std::list<EdgeWPtr>::iterator Edge::getVertexDstListIt() const {
    return this->vertex_dst_list_it_;
}

void Edge::setVertexDstListIt(std::list<EdgeWPtr>::iterator list_it) {
    this->vertex_dst_list_it_ = list_it;
}

FacetSPtr Edge::getFacetL() const {
    // DEBUG_WPTR(this->facet_l_);
    if (this->facet_l_.expired())
        return FacetSPtr();
    else
        return FacetSPtr(this->facet_l_);
}

void Edge::setFacetL(FacetSPtr facet) {
    this->facet_l_ = facet;
}

std::list<EdgeSPtr>::iterator Edge::getFacetLListIt() const {
    return this->facet_l_list_it_;
}

void Edge::setFacetLListIt(std::list<EdgeSPtr>::iterator list_it) {
    this->facet_l_list_it_ = list_it;
}

FacetSPtr Edge::getFacetR() const {
    // DEBUG_WPTR(this->facet_r_);
    if (this->facet_r_.expired())
        return FacetSPtr();
    else
        return FacetSPtr(this->facet_r_);
}

void Edge::setFacetR(FacetSPtr facet) {
    this->facet_r_ = facet;
}

std::list<EdgeSPtr>::iterator Edge::getFacetRListIt() const {
    return this->facet_r_list_it_;
}

void Edge::setFacetRListIt(std::list<EdgeSPtr>::iterator list_it) {
    this->facet_r_list_it_ = list_it;
}

PolyhedronSPtr Edge::getPolyhedron() const {
    // DEBUG_WPTR(this->polyhedron_);
    if (this->polyhedron_.expired())
        return PolyhedronSPtr();
    else
        return PolyhedronSPtr(this->polyhedron_);
}

void Edge::setPolyhedron(PolyhedronSPtr polyhedron) {
    this->polyhedron_ = polyhedron;
}

std::list<EdgeSPtr>::iterator Edge::getPolyhedronListIt() const {
    return this->polyhedron_list_it_;
}

void Edge::setPolyhedronListIt(std::list<EdgeSPtr>::iterator list_it) {
    this->polyhedron_list_it_ = list_it;
}

EdgeDataSPtr Edge::getData() const {
    DEBUG_SPTR(this->data_);
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

Segment3SPtr Edge::segment() const {
    return KernelFactory::createSegment3(
            vertex_src_->getPoint(), vertex_dst_->getPoint());
}

Line3SPtr Edge::line() const {
    return KernelFactory::createLine3(
            vertex_src_->getPoint(), vertex_dst_->getPoint());
}

FacetSPtr Edge::other(FacetSPtr facet) const {
    FacetSPtr result = FacetSPtr();
    if (facet == facet_l_.lock()) {
        result = facet_r_.lock();
    } else if (facet == facet_r_.lock()) {
        result = facet_l_.lock();
    }
    return result;
}

VertexSPtr Edge::src(FacetSPtr facet_l) const {
    VertexSPtr result = VertexSPtr();
    if (facet_l == facet_l_.lock()) {
        result = vertex_src_;
    } else if (facet_l == facet_r_.lock()) {
        result = vertex_dst_;
    }
    return result;
}

VertexSPtr Edge::dst(FacetSPtr facet_l) const {
    VertexSPtr result = VertexSPtr();
    if (facet_l == facet_l_.lock()) {
        result = vertex_dst_;
    } else if (facet_l == facet_r_.lock()) {
        result = vertex_src_;
    }
    return result;
}

FacetSPtr Edge::left(VertexSPtr vertex_src) const {
    FacetSPtr result = FacetSPtr();
    if (vertex_src == vertex_src_) {
        if (!facet_l_.expired()) {
            result = FacetSPtr(facet_l_);
        }
    } else if (vertex_src == vertex_dst_) {
        if (!facet_r_.expired()) {
            result = FacetSPtr(facet_r_);
        }
    }
    return result;
}

FacetSPtr Edge::right(VertexSPtr vertex_src) const {
    FacetSPtr result = FacetSPtr();
    if (vertex_src == vertex_src_) {
        if (!facet_r_.expired()) {
            result = FacetSPtr(facet_r_);
        }
    } else if (vertex_src == vertex_dst_) {
        if (!facet_l_.expired()) {
            result = FacetSPtr(facet_l_);
        }
    }
    return result;
}


EdgeSPtr Edge::next(FacetSPtr facet) const {
    EdgeSPtr result = EdgeSPtr();
    std::list<EdgeSPtr>::const_iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e;
        if (edge == shared_from_this()) {
            result = edge;
            break;
        }
        it_e++;
    }
    if (it_e != facet->edges().end()) {
        VertexSPtr vertex_dst = this->dst(facet);
        std::list<EdgeSPtr>::const_iterator it_e_begin = it_e++;
        if (it_e == facet->edges().end()) {
            it_e = facet->edges().begin();
        }
        while (it_e != it_e_begin) {
            EdgeSPtr edge = *it_e++;
            if (it_e == facet->edges().end()) {
                it_e = facet->edges().begin();
            }
            if (vertex_dst == edge->src(facet)) {
                result = edge;
                break;
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

EdgeSPtr Edge::prev(FacetSPtr facet) const {
    EdgeSPtr result = EdgeSPtr();
    std::list<EdgeSPtr>::const_reverse_iterator it_e = facet->edges().rbegin();
    while (it_e != facet->edges().rend()) {
        EdgeSPtr edge = *it_e;
        if (edge == shared_from_this()) {
            result = edge;
            break;
        }
        it_e++;
    }
    if (it_e != facet->edges().rend()) {
        VertexSPtr vertex_src = this->src(facet);
        std::list<EdgeSPtr>::const_reverse_iterator it_e_begin = it_e++;
        if (it_e == facet->edges().rend()) {
            it_e = facet->edges().rbegin();
        }
        while (it_e != it_e_begin) {
            EdgeSPtr edge = *it_e++;
            if (it_e == facet->edges().rend()) {
                it_e = facet->edges().rbegin();
            }
            if (vertex_src == edge->dst(facet)) {
                result = edge;
                break;
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

EdgeSPtr Edge::next(VertexSPtr vertex) const {
    EdgeSPtr result = EdgeSPtr();
    FacetSPtr facet;
    if (vertex == vertex_src_) {
        facet = FacetSPtr(facet_l_);
    } else if (vertex == vertex_dst_) {
        facet = FacetSPtr(facet_r_);
    }
    if (facet) {
        std::list<EdgeSPtr> edges_possible;
        std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
        while (it_e != vertex->edges().end()) {
            EdgeWPtr edge_wptr = *it_e++;
            if (!edge_wptr.expired()) {
                EdgeSPtr edge(edge_wptr);
                if (edge.get() == this) {
                    continue;
                }
                if (edge->dst(facet) == vertex) {
                    edges_possible.push_back(edge);
                }
            }
        }
        if (edges_possible.size() == 1) {
            result = edges_possible.front();
        } else {
            double angle_min = 2*M_PI;
            std::list<EdgeSPtr>::iterator it_e = edges_possible.begin();
            while (it_e != edges_possible.end()) {
                EdgeSPtr edge = *it_e++;
                double angle = angleTo(edge);
                 if (angle == angle_min) {
                    DEBUG_VAL("Warning: Not able to distinguish possible next edges.");
                }
                if (angle <= angle_min) {
                    result = edge;
                    angle_min = angle;
                }
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

EdgeSPtr Edge::prev(VertexSPtr vertex) const {
    EdgeSPtr result = EdgeSPtr();
    FacetSPtr facet;
    if (vertex == vertex_src_) {
        facet = FacetSPtr(facet_r_);
    } else if (vertex == vertex_dst_) {
        facet = FacetSPtr(facet_l_);
    }
    if (facet) {
        std::list<EdgeSPtr> edges_possible;
        std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
        while (it_e != vertex->edges().end()) {
            EdgeWPtr edge_wptr = *it_e++;
            if (!edge_wptr.expired()) {
                EdgeSPtr edge(edge_wptr);
                if (edge.get() == this) {
                    continue;
                }
                if (edge->src(facet) == vertex) {
                    edges_possible.push_back(edge);
                }
            }
        }
        if (edges_possible.size() == 1) {
            result = edges_possible.front();
        } else {
            double angle_max = 0.0;
            std::list<EdgeSPtr>::iterator it_e = edges_possible.begin();
            while (it_e != edges_possible.end()) {
                EdgeSPtr edge = *it_e++;
                double angle = angleTo(edge);
                if (angle == angle_max) {
                    DEBUG_VAL("Warning: Not able to distinguish possible next edges.");
                }
                if (angle >= angle_max) {
                    result = edge;
                    angle_max = angle;
                }
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

void Edge::invert() {
    VertexSPtr vertex_tmp_ = vertex_src_;
    std::list<EdgeWPtr>::iterator vertex_tmp_list_it_ = vertex_src_list_it_;
    FacetWPtr facet_tmp_ = facet_l_;
    std::list<EdgeSPtr>::iterator facet_tmp_list_it_ = facet_l_list_it_;

    vertex_src_ = vertex_dst_;
    vertex_src_list_it_ = vertex_dst_list_it_;
    facet_l_ = facet_r_;
    facet_l_list_it_ = facet_r_list_it_;

    vertex_dst_ = vertex_tmp_;
    vertex_dst_list_it_ = vertex_tmp_list_it_;
    facet_r_ = facet_tmp_;
    facet_r_list_it_ = facet_tmp_list_it_;
}

EdgeSPtr Edge::split(VertexSPtr middle) {
    EdgeSPtr result = Edge::create(middle, vertex_dst_);
    if (!facet_l_.expired()) {
        FacetSPtr facet_l = FacetSPtr(facet_l_);
        result->setFacetL(facet_l);
        facet_l->addEdge(result);
    }
    if (!facet_r_.expired()) {
        FacetSPtr facet_r = FacetSPtr(facet_r_);
        result->setFacetR(facet_r);
        facet_r->addEdge(result);
    }
    vertex_dst_->removeEdge(shared_from_this());
    vertex_dst_ = middle;
    middle->addEdge(shared_from_this());
    if (!polyhedron_.expired()) {
        PolyhedronSPtr polyhedron = PolyhedronSPtr(polyhedron_);
        polyhedron->addEdge(result);
    }
    return result;
}

void Edge::replaceVertexSrc(VertexSPtr vertex_src) {
    vertex_src_->removeEdge(shared_from_this());
    vertex_src_ = vertex_src;
    vertex_src->addEdge(shared_from_this());
}

void Edge::replaceVertexDst(VertexSPtr vertex_dst) {
    vertex_dst_->removeEdge(shared_from_this());
    vertex_dst_ = vertex_dst;
    vertex_dst->addEdge(shared_from_this());
}

void Edge::replaceFacetL(FacetSPtr facet_l) {
    if (!facet_l_.expired()) {
        FacetSPtr facet = FacetSPtr(facet_l_);
        facet->removeEdge(shared_from_this());
    }
    facet_l_ = facet_l;
    facet_l->addEdge(shared_from_this());
}

void Edge::replaceFacetR(FacetSPtr facet_r) {
    if (!facet_r_.expired()) {
        FacetSPtr facet = FacetSPtr(facet_r_);
        facet->removeEdge(shared_from_this());
    }
    facet_r_ = facet_r;
    facet_r->addEdge(shared_from_this());
}

bool Edge::hasSameFacets(EdgeSPtr edge) {
    bool result = (
        (facet_r_ == edge->facet_r_ && facet_l_ == edge->facet_l_) ||
        (facet_r_ == edge->facet_l_ && facet_l_ == edge->facet_r_));
    return result;
}

int Edge::getID() const {
    return this->id_;
}

void Edge::setID(int id) {
    this->id_ = id;
}

double Edge::angle() const {
    double result = 0.0;
    if (!facet_l_.expired() && !facet_r_.expired()) {
        FacetSPtr facet_l = facet_l_.lock();
        FacetSPtr facet_r = facet_r_.lock();
        Vector3SPtr v1 = KernelFactory::createVector3(facet_l->plane());
        Vector3SPtr v2 = KernelFactory::createVector3(facet_r->plane());
#ifdef USE_CGAL
        result = acos(((*v1) * (*v2)) /
                CGAL::sqrt(v1->squared_length() * v2->squared_length()));
#else
        result = acos(((*v1) * (*v2)) /
                sqrt(v1->squared_length() * v2->squared_length()));
#endif
        result = M_PI - result;
        if (isReflex()) {
            result = 2.0*M_PI - result;
        }
    } else {
        DEBUG_VAL("Warning: Not able to determine angle.");
        DEBUG_VAL(toString());
    }
    return result;
}

bool Edge::isReflex() const {
    bool result = false;
    if (vertex_src_->getPoint() != vertex_dst_->getPoint() &&
            !facet_l_.expired() && !facet_r_.expired()) {
        FacetSPtr facet_l = FacetSPtr(facet_l_);
        FacetSPtr facet_r = FacetSPtr(facet_r_);
        Plane3SPtr plane_l = facet_l->plane();
        Plane3SPtr plane_r = facet_r->plane();
        Vector3SPtr dir = KernelFactory::createVector3(line());
        Vector3SPtr normal_l = KernelFactory::createVector3(plane_l);
        Point3SPtr p_src = vertex_src_->getPoint();
#ifdef USE_CGAL
        Point3 p = (*p_src) + CGAL::cross_product(*normal_l, *dir);
        if (plane_r->oriented_side(p) == CGAL::ON_POSITIVE_SIDE) {
            result = true;
        }
#else
        Point3 p = (*p_src) + normal_l->cross(*dir);
        if (plane_r->side(p) > 0) {
            result = true;
        }
#endif
    } else {
        DEBUG_VAL("Warning: Not able to determine if edge is reflex.");
        DEBUG_VAL(toString());
    }
    return result;
}


double Edge::angleTo(EdgeSPtr edge) const {
    double result = 0.0;
    VertexSPtr vertex;
    if (vertex_src_ == edge->getVertexSrc() ||
            vertex_src_ == edge->getVertexDst()) {
        vertex = vertex_src_;
    } else if (vertex_dst_ == edge->getVertexSrc() ||
            vertex_dst_ == edge->getVertexDst()) {
        vertex = vertex_dst_;
    }
    FacetSPtr facet;
    FacetSPtr facet_l = getFacetL();
    FacetSPtr facet_r = getFacetR();
    if (facet_l == edge->getFacetL() ||
            facet_l == edge->getFacetR()) {
        facet = facet_l;
    } else if (facet_r == edge->getFacetL() ||
            facet_r == edge->getFacetR()) {
        facet = facet_r;
    }
    if (vertex && facet) {
        Vector3SPtr normal = KernelFactory::createVector3(facet->plane());
        Vector3SPtr dir_self;
        if (vertex == vertex_src_) {
            dir_self = KernelFactory::createVector3(
                    *(vertex_dst_->getPoint()) - *(vertex_src_->getPoint()));
        } else {
            dir_self = KernelFactory::createVector3(
                    *(vertex_src_->getPoint()) - *(vertex_dst_->getPoint()));
        }
        Vector3SPtr dir_other;
        if (vertex == edge->getVertexSrc()) {
            dir_other = KernelFactory::createVector3(
                    *(edge->getVertexDst()->getPoint()) - *(edge->getVertexSrc()->getPoint()));
        } else {
            dir_other = KernelFactory::createVector3(
                    *(edge->getVertexSrc()->getPoint()) - *(edge->getVertexDst()->getPoint()));
        }
        double angle = 0.0;
#ifdef USE_CGAL
        angle = acos(((*dir_self) * (*dir_other)) /
                CGAL::sqrt(dir_self->squared_length() * dir_other->squared_length()));
#else
        angle = acos(((*dir_self) * (*dir_other)) /
                sqrt(dir_self->squared_length() * dir_other->squared_length()));
#endif
        if ((facet_l == edge->getFacetR() && facet_r == edge->getFacetL()) ||
                (facet_l == edge->getFacetL() && facet_r == edge->getFacetR())) {
            if (angle < M_PI/2.0) {
                result = 0.0;
            } else {
                result = M_PI;
            }
        } else {
            result = angle;
            double angle_normal = 0.0;
#ifdef USE_CGAL
            Vector3 crossprod = CGAL::cross_product(*dir_self, *dir_other);
            angle_normal = acos(((*normal) * crossprod) /
                    CGAL::sqrt(normal->squared_length() * crossprod.squared_length()));
#else
            Vector3 crossprod = dir_self->cross(*dir_other);
            angle_normal = acos(((*normal) * crossprod) /
                    sqrt(normal->squared_length() * crossprod.squared_length()));
#endif
            if (angle_normal > M_PI/2.0) {
                result += M_PI;
            }
        }
    } else {
        DEBUG_PRINT("Warning: Edges do not have a shared Vertex and a shared Facet.");
    }
    return result;
}


std::string Edge::toString() const {
    std::string result("Edge(");
    if (id_ != -1) {
        result += "id=" + util::StringFactory::fromInteger(id_) + ", ";
    } else {
        result += util::StringFactory::fromPointer(this) + ", ";
    }
    result += "src=" + vertex_src_->toString() + ", ";
    result += "dst=" + vertex_dst_->toString();
    if (!facet_l_.expired()) {
        result += ", l=";
        if (getFacetL()->getID() != -1) {
            result += util::StringFactory::fromInteger(getFacetL()->getID());
        } else {
            result += util::StringFactory::fromPointer(getFacetL().get());
        }
    }
    if (!facet_r_.expired()) {
        result += ", r=";
        if (getFacetR()->getID() != -1) {
            result += util::StringFactory::fromInteger(getFacetR()->getID());
        } else {
            result += util::StringFactory::fromPointer(getFacetR().get());
        }
    }
    result += ")";
    return result;
}

} }
