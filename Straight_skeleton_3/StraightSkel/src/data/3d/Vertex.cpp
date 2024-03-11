/**
 * @file   data/3d/Vertex.cpp
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#include "data/3d/Vertex.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "util/StringFactory.h"
#include <algorithm>

namespace data { namespace _3d {

Vertex::Vertex(Point3SPtr point) {
    this->point_ = point;
    this->id_ = -1;
}

Vertex::~Vertex() {
    facets_.clear();
    edges_.clear();
}

Vertex::Vertex(const Vertex& vertex) {
    point_ = vertex.point_;
    this->id_ = -1;
}

VertexSPtr Vertex::create(Point3SPtr point) {
    VertexSPtr result = VertexSPtr(new Vertex(point));
    return result;
}

VertexSPtr Vertex::clone() const {
    VertexSPtr result = VertexSPtr(new Vertex(*this));
    return result;
}

Point3SPtr Vertex::getPoint() const {
    DEBUG_SPTR(this->point_);
    return this->point_;
}

void Vertex::setPoint(Point3SPtr point) {
    this->point_ = point;
}

void Vertex::addEdge(EdgeSPtr edge) {
    EdgeWPtr edge_wptr(edge);
    std::list<EdgeWPtr>::iterator it = edges_.insert(edges_.end(), edge_wptr);
    VertexSPtr vertex_src = edge->getVertexSrc();
    VertexSPtr vertex_dst = edge->getVertexDst();
    if (vertex_src == shared_from_this() && vertex_dst == shared_from_this()) {
        std::list<EdgeWPtr>::iterator it_e =
                std::find(edges_.begin(), edges_.end(), edge_wptr);
        if (it_e == edge->getVertexSrcListIt()) {
            edge->setVertexDstListIt(it);
        } else {
            edge->setVertexSrcListIt(it);
        }
    } else if (vertex_src == shared_from_this()) {
        edge->setVertexSrcListIt(it);
    } else if (vertex_dst == shared_from_this()) {
        edge->setVertexDstListIt(it);
    }
}

bool Vertex::removeEdge(EdgeSPtr edge) {
    bool result = false;
    if (edge->getVertexSrc() == shared_from_this()) {
        edges_.erase(edge->getVertexSrcListIt());
        edge->setVertexSrcListIt(std::list<EdgeWPtr>::iterator());
        result = true;
    } else if (edge->getVertexDst() == shared_from_this()) {
        edges_.erase(edge->getVertexDstListIt());
        edge->setVertexDstListIt(std::list<EdgeWPtr>::iterator());
        result = true;
    }
    return result;
}

EdgeSPtr Vertex::firstEdge() const {
    EdgeSPtr result;
    std::list<EdgeWPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            result = EdgeSPtr(edge_wptr);
            break;
        }
    }
    DEBUG_SPTR(result);
    return result;
}

EdgeSPtr Vertex::getEdge(unsigned int i) {
    EdgeSPtr edge = firstEdge();
    for (unsigned int j = 0; j < i; j++) {
        edge = edge->next(shared_from_this());
    }
    return edge;
}

EdgeSPtr Vertex::findEdge(VertexSPtr dst) const {
    EdgeSPtr result = EdgeSPtr();
    std::list<EdgeWPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge = EdgeSPtr(edge_wptr);
            if (edge->getVertexSrc().get() == this &&
                    edge->getVertexDst() == dst) {
                result = edge;
                break;
            }
            if (edge->getVertexDst().get() == this &&
                    edge->getVertexSrc() == dst) {
                result = edge;
                break;
            }
        }
    }
    return result;
}

EdgeSPtr Vertex::findEdge(FacetSPtr facet) const {
    EdgeSPtr result = EdgeSPtr();
    std::list<EdgeWPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge = EdgeSPtr(edge_wptr);
            if (edge->src(facet).get() == this) {
                result = edge;
                break;
            }
        }
    }
    return result;
}

void Vertex::addFacet(FacetSPtr facet) {
    facets_.insert(facets_.end(), FacetWPtr(facet));
}

bool Vertex::removeFacet(FacetSPtr facet) {
    bool result = false;
    std::list<FacetWPtr>::iterator it = facets_.begin();
    while (it != facets_.end()) {
        std::list<FacetWPtr>::iterator it_current = it;
        FacetWPtr facet_wptr = *it++;
        if (!facet_wptr.expired()) {
            if (facet_wptr.lock() == facet) {
                facets_.erase(it_current);
                result = true;
                break;
            }
        }
    }
    return result;
}

FacetSPtr Vertex::firstFacet() const {
    FacetSPtr result;
    std::list<FacetWPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (!facet_wptr.expired()) {
            result = FacetSPtr(facet_wptr);
            break;
        }
    }
    DEBUG_SPTR(result);
    return result;
}

FacetSPtr Vertex::getFacet(unsigned int i) {
    FacetSPtr facet = firstFacet();
    for (unsigned int j = 0; j < i; j++) {
        facet = facet->next(shared_from_this());
    }
    return facet;
}

bool Vertex::containsEdge(EdgeSPtr edge) const {
    EdgeWPtr edge_wptr = EdgeWPtr(edge);
    bool result = (edges_.end() !=
            std::find(edges_.begin(), edges_.end(), edge_wptr));
    return result;
}

bool Vertex::containsFacet(FacetSPtr facet) const {
    FacetWPtr facet_wptr = FacetWPtr(facet);
    bool result = (facets_.end() !=
            std::find(facets_.begin(), facets_.end(), facet_wptr));
    return result;
}

void Vertex::sortEdges() {
    std::list<EdgeWPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        std::list<EdgeWPtr>::iterator it_current = it_e;
        EdgeWPtr edge_wptr = *it_e++;
        if (edge_wptr.expired()) {
            edges_.erase(it_current);
        }
    }
    std::list<EdgeSPtr> tmp;
    std::list<EdgeSPtr>::iterator it_e_tmp = tmp.begin();
    if (edges_.size() > 0) {
        EdgeSPtr first = EdgeSPtr();
        EdgeSPtr edge = EdgeSPtr(edges_.front());
        while (edge != first) {
            if (!first) {
                first = edge;
            }
            tmp.insert(tmp.end(), edge);
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

void Vertex::sortFacets() {
    std::list<FacetWPtr>::iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        std::list<FacetWPtr>::iterator it_current = it_f;
        FacetWPtr facet_wptr = *it_f++;
        if (facet_wptr.expired()) {
            facets_.erase(it_current);
        }
    }
    std::list<FacetSPtr> tmp;
    std::list<FacetSPtr>::iterator it_f_tmp = tmp.begin();
    if (facets_.size() > 0) {
        FacetSPtr facet = FacetSPtr(facets_.front());
        EdgeSPtr edge_first = findEdge(facet);
        EdgeSPtr edge;
        while (edge != edge_first) {
            if (!edge) {
                edge = edge_first;
            }
            tmp.insert(tmp.end(), facet);
            edge = edge->next(shared_from_this());
            facet = edge->other(facet);
        }
        while (it_f_tmp != tmp.end()) {
            FacetSPtr facet = *it_f_tmp++;
            removeFacet(facet);
        }
    }
    facets_.clear();
    it_f_tmp = tmp.begin();
    while (it_f_tmp != tmp.end()) {
        FacetSPtr facet = *it_f_tmp++;
        addFacet(facet);
    }
}

void Vertex::sort() {
    std::list<EdgeWPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        std::list<EdgeWPtr>::iterator it_current = it_e;
        EdgeWPtr edge_wptr = *it_e++;
        if (edge_wptr.expired()) {
            edges_.erase(it_current);
        }
    }
    std::list<FacetWPtr>::iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        std::list<FacetWPtr>::iterator it_current = it_f;
        FacetWPtr facet_wptr = *it_f++;
        if (facet_wptr.expired()) {
            facets_.erase(it_current);
        }
    }
    std::list<EdgeSPtr> edges_tmp;
    std::list<EdgeSPtr>::iterator it_e_tmp = edges_tmp.begin();
    std::list<FacetSPtr> facets_tmp;
    std::list<FacetSPtr>::iterator it_f_tmp = facets_tmp.begin();
    if (edges_.size() > 0) {
        EdgeSPtr edge_first = EdgeSPtr(edges_.front());
        EdgeSPtr edge;
        FacetSPtr facet = edge_first->getFacetL();
        if (edge_first->getVertexDst() == shared_from_this()) {
            facet = edge_first->getFacetR();
        }
        while (edge != edge_first) {
            if (!edge) {
                edge = edge_first;
            }
            edges_tmp.insert(edges_tmp.end(), edge);
            facets_tmp.insert(facets_tmp.end(), facet);
            edge = edge->next(shared_from_this());
            facet = edge->other(facet);
        }
        while (it_e_tmp != edges_tmp.end()) {
            EdgeSPtr edge = *it_e_tmp++;
            removeEdge(edge);
        }
        while (it_f_tmp != facets_tmp.end()) {
            FacetSPtr facet = *it_f_tmp++;
            removeFacet(facet);
        }
    }
    edges_.clear();
    it_e_tmp = edges_tmp.begin();
    while (it_e_tmp != edges_tmp.end()) {
        EdgeSPtr edge = *it_e_tmp++;
        addEdge(edge);
    }
    facets_.clear();
    it_f_tmp = facets_tmp.begin();
    while (it_f_tmp != facets_tmp.end()) {
        FacetSPtr facet = *it_f_tmp++;
        addFacet(facet);
    }
}

PolyhedronSPtr Vertex::getPolyhedron() const {
    // DEBUG_WPTR(this->polyhedron_);
    if (this->polyhedron_.expired())
        return PolyhedronSPtr();
    else
        return PolyhedronSPtr(this->polyhedron_);
}

void Vertex::setPolyhedron(PolyhedronSPtr polyhedron) {
    this->polyhedron_ = polyhedron;
}

std::list<VertexSPtr>::iterator Vertex::getPolyhedronListIt() const {
    return this->polyhedron_list_it_;
}

void Vertex::setPolyhedronListIt(std::list<VertexSPtr>::iterator list_it) {
    this->polyhedron_list_it_ = list_it;
}

VertexDataSPtr Vertex::getData() const {
    DEBUG_SPTR(this->data_);
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

std::list<EdgeWPtr>& Vertex::edges() {
    return this->edges_;
}

std::list<FacetWPtr>& Vertex::facets() {
    return this->facets_;
}

VertexSPtr Vertex::next(FacetSPtr facet) const {
    VertexSPtr result = VertexSPtr();
    std::list<VertexSPtr>::const_iterator it_v = facet->vertices().begin();
    while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v;
        if (vertex == shared_from_this()) {
            if (facet->vertices().size() == 1) {
                result = vertex;
            }
            break;
        }
        it_v++;
    }
    if (it_v != facet->vertices().end()) {
        std::list<VertexSPtr>::const_iterator it_v_begin = it_v++;
        if (it_v == facet->vertices().end()) {
            it_v = facet->vertices().begin();
        }
        while (it_v != it_v_begin) {
            VertexSPtr vertex = *it_v++;
            if (it_v == facet->vertices().end()) {
                it_v = facet->vertices().begin();
            }
            EdgeSPtr edge = findEdge(vertex);
            if (edge) {
                if (edge->dst(facet) == vertex) {
                    result = vertex;
                    break;
                }
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

VertexSPtr Vertex::prev(FacetSPtr facet) const {
    VertexSPtr result = VertexSPtr();
    std::list<VertexSPtr>::const_reverse_iterator it_v = facet->vertices().rbegin();
    while (it_v != facet->vertices().rend()) {
        VertexSPtr vertex = *it_v;
        if (vertex == shared_from_this()) {
            if (facet->vertices().size() == 1) {
                result = vertex;
            }
            break;
        }
        it_v++;
    }
    if (it_v != facet->vertices().rend()) {
        std::list<VertexSPtr>::const_reverse_iterator it_v_begin = it_v++;
        if (it_v == facet->vertices().rend()) {
            it_v = facet->vertices().rbegin();
        }
        while (it_v != it_v_begin) {
            VertexSPtr vertex = *it_v++;
            if (it_v == facet->vertices().rend()) {
                it_v = facet->vertices().rbegin();
            }
            EdgeSPtr edge = findEdge(vertex);
            if (edge) {
                if (edge->src(facet) == vertex) {
                    result = vertex;
                    break;
                }
            }
        }
    }
    DEBUG_SPTR(result);
    return result;
}

VertexSPtr Vertex::split(FacetSPtr facet_left, FacetSPtr facet_right) {
    VertexSPtr result = VertexSPtr();
    if (!containsFacet(facet_left) || !containsFacet(facet_right)) {
        return result;
    }
//    if (facet_left->findEdge(facet_right)) {
//        return result;
//    }
    result = Vertex::create(getPoint());

    // 1. select edges and facets for result
    std::list<FacetSPtr> facets;
    std::list<EdgeSPtr> edges;
    FacetSPtr poly_curr = facet_right;
    EdgeSPtr edge_curr = findEdge(poly_curr);
    while (poly_curr != facet_left) {
        facets.insert(facets.end(), poly_curr);
        edge_curr = edge_curr->next(shared_from_this());
        edges.insert(edges.end(), edge_curr);
        poly_curr = edge_curr->other(poly_curr);
    }
    facets.insert(facets.end(), facet_left);

    // 2. split vertex
    std::list<EdgeSPtr>::iterator it_e = edges.begin();
    while (it_e != edges.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getVertexSrc() == shared_from_this()) {
            edge->replaceVertexSrc(result);
        } else if (edge->getVertexDst() == shared_from_this()) {
            edge->replaceVertexDst(result);
        }
    }
    std::list<FacetSPtr>::iterator it_f = facets.begin();
    while (it_f != facets.end()) {
        FacetSPtr facet = *it_f++;
        if (facet != facet_left && facet != facet_right) {
            facet->removeVertex(shared_from_this());
        }
        facet->addVertex(result);
    }

    // 3. insert connecting edge
    EdgeSPtr edge = Edge::create(shared_from_this(), result);
    edge->setFacetL(facet_left);
    edge->setFacetR(facet_right);
    facet_left->addEdge(edge);
    facet_right->addEdge(edge);

    if (!polyhedron_.expired()) {
        PolyhedronSPtr polyhedron = PolyhedronSPtr(polyhedron_);
        polyhedron->addVertex(result);
        polyhedron->addEdge(edge);
    }
    return result;
}

double Vertex::getX() const {
#ifdef USE_CGAL
    return this->point_->x();
#else
    return this->point_->getX();
#endif
}

double Vertex::getY() const {
#ifdef USE_CGAL
    return this->point_->y();
#else
    return this->point_->getY();
#endif
}

double Vertex::getZ() const {
#ifdef USE_CGAL
    return this->point_->z();
#else
    return this->point_->getZ();
#endif
}

int Vertex::getID() const {
    return this->id_;
}

void Vertex::setID(int id) {
    this->id_ = id;
}

unsigned int Vertex::degree() const {
    unsigned int result = 0;
    std::list<EdgeWPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            result++;
        }
    }
    return result;
}

bool Vertex::isReflex() const {
    if (degree() == 0) {
        return false;
    }
    bool result = true;
    std::list<EdgeWPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge = EdgeSPtr(edge_wptr);
            if (!edge->isReflex()) {
                result = false;
            }
        }
    }
    return result;
}

bool Vertex::isConvex() const {
    if (degree() == 0) {
        return false;
    }
    bool result = true;
    std::list<EdgeWPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeWPtr edge_wptr = *it_e++;
        if (!edge_wptr.expired()) {
            EdgeSPtr edge = EdgeSPtr(edge_wptr);
            if (edge->isReflex()) {
                result = false;
            }
        }
    }
    return result;
}

std::string Vertex::toString() const {
    std::string result("Vertex(");
    if (id_ != -1) {
        result += "id=" + util::StringFactory::fromInteger(id_) + ", ";
    } else {
        result += util::StringFactory::fromPointer(this) + ", ";
    }
    result += "<" + util::StringFactory::fromDouble(getX()) + ", ";
    result += util::StringFactory::fromDouble(getY()) + ", ";
    result += util::StringFactory::fromDouble(getZ()) + ">";
//    if (edges_.size() > 0) {
//        result += ", Edges:" + util::StringFactory::fromInteger(edges_.size());
//    }
//    if (facets_.size() > 0) {
//        result += ", Facets:" + util::StringFactory::fromInteger(facets_.size());
//    }
    result += ")";
    return result;
}

} }
