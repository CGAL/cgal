/**
 * @file   data/3d/SphericalPolygon.cpp
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#include "data/3d/SphericalPolygon.h"

#include "debug.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/KernelFactory.h"
#include "util/StringFactory.h"
#include <map>
#include <sstream>

namespace data { namespace _3d {

SphericalPolygon::SphericalPolygon(Sphere3SPtr sphere) {
    this->sphere_ = sphere;
}

SphericalPolygon::~SphericalPolygon() {
    sphere_.reset();
    edges_.clear();
    vertices_.clear();
}


SphericalPolygonSPtr SphericalPolygon::create(Sphere3SPtr sphere) {
    SphericalPolygonSPtr result = SphericalPolygonSPtr(new SphericalPolygon(sphere));
    return result;
}

SphericalPolygonSPtr SphericalPolygon::create(Sphere3SPtr sphere,
        unsigned int num_edges, CircularEdgeSPtr edges[]) {
    SphericalPolygonSPtr result = SphericalPolygonSPtr(new SphericalPolygon(sphere));
    for (unsigned int i = 0; i < num_edges; i++) {
        result->addEdge(edges[i]);
    }
    return result;
}


SphericalPolygonSPtr SphericalPolygon::clone() const {
    std::map<CircularVertexSPtr, CircularVertexSPtr> vertices_c;
    SphericalPolygonSPtr result = SphericalPolygon::create(sphere_);
    std::list<CircularVertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        CircularVertexSPtr vertex = *it_v++;
        CircularVertexSPtr vertex_c = vertex->clone();
        result->addVertex(vertex_c);
        vertices_c[vertex] = vertex_c;
    }
    std::list<CircularEdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        CircularEdgeSPtr edge = *it_e++;
        CircularVertexSPtr src = vertices_c[edge->getVertexSrc()];
        CircularVertexSPtr dst = vertices_c[edge->getVertexDst()];
        CircularEdgeSPtr edge_c = CircularEdge::create(src, dst);
        result->addEdge(edge_c);
    }
    return result;
}


Sphere3SPtr SphericalPolygon::getSphere() const {
    return sphere_;
}

void SphericalPolygon::addVertex(CircularVertexSPtr vertex) {
    std::list<CircularVertexSPtr>::iterator it = vertices_.insert(vertices_.end(), vertex);
    vertex->setPolygon(shared_from_this());
    vertex->setListIt(it);
}

bool SphericalPolygon::removeVertex(CircularVertexSPtr vertex) {
    bool result = false;
    if (vertex->getPolygon() == shared_from_this()) {
        vertices_.erase(vertex->getListIt());
        vertex->setPolygon(SphericalPolygonSPtr());
        vertex->setListIt(std::list<CircularVertexSPtr>::iterator());
        CircularEdgeSPtr edge = vertex->getEdgeIn();
        if (edge) {
            this->removeEdge(edge);
            vertex->setEdgeIn(CircularEdgeSPtr());
        }
        edge = vertex->getEdgeOut();
        if (edge) {
            this->removeEdge(edge);
            vertex->setEdgeOut(CircularEdgeSPtr());
        }
        result = true;
    }
    return result;
}

void SphericalPolygon::addEdge(CircularEdgeSPtr edge) {
    std::list<CircularEdgeSPtr>::iterator it = edges_.insert(edges_.end(), edge);
    edge->setPolygon(shared_from_this());
    edge->setListIt(it);
    CircularVertexSPtr vertex = edge->getVertexSrc();
    if (vertex->getPolygon() != shared_from_this()) {
        this->addVertex(vertex);
    }
    vertex = edge->getVertexDst();
    if (vertex->getPolygon() != shared_from_this()) {
        this->addVertex(vertex);
    }
}

bool SphericalPolygon::removeEdge(CircularEdgeSPtr edge) {
    bool result = false;
    if (edge->getPolygon() == shared_from_this()) {
        edges_.erase(edge->getListIt());
        edge->setPolygon(SphericalPolygonSPtr());
        edge->setListIt(std::list<CircularEdgeSPtr>::iterator());
        CircularVertexSPtr vertex = edge->getVertexSrc();
        if (vertex->getEdgeOut() == edge) {
            vertex->setEdgeOut(CircularEdgeSPtr());
        }
        vertex = edge->getVertexDst();
        if (vertex->getEdgeIn() == edge) {
            vertex->setEdgeIn(CircularEdgeSPtr());
        }
        result = true;
    }
    return result;
}

void SphericalPolygon::clear() {
    std::list<CircularEdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        CircularEdgeSPtr edge = *it_e++;
        this->removeEdge(edge);
    }
    std::list<CircularVertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        CircularVertexSPtr vertex = *it_v++;
        this->removeVertex(vertex);
    }
}


SharedMutex& SphericalPolygon::mutex() {
    return this->mutex_;
}

std::list<CircularVertexSPtr>& SphericalPolygon::vertices() {
    return this->vertices_;
}

std::list<CircularEdgeSPtr>& SphericalPolygon::edges() {
    return this->edges_;
}

bool SphericalPolygon::isConsistent() const  {
    bool result = true;
    std::list<CircularVertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        CircularVertexSPtr vertex = *it_v++;
        if (vertex->getPolygon() != shared_from_this()) {
            DEBUG_VAR(vertex->toString());
            result = false;
        }
        if (vertex->getEdgeOut()) {
            if (vertex->getEdgeOut()->getVertexSrc() != vertex) {
                DEBUG_VAR(vertex->toString());
                DEBUG_VAR(vertex->getEdgeOut()->toString());
                result = false;
            }
        }
        if (vertex->getEdgeIn()) {
            if (vertex->getEdgeIn()->getVertexDst() != vertex) {
                DEBUG_VAR(vertex->toString());
                DEBUG_VAR(vertex->getEdgeOut()->toString());
                result = false;
            }
        }
    }
    std::list<CircularEdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        CircularEdgeSPtr edge = *it_e++;
        if (edge->getPolygon() != shared_from_this()) {
            DEBUG_VAR(edge->toString());
            result = false;
        }
        if (edge->getVertexSrc()) {
            if (edge->getVertexSrc()->getEdgeOut() != edge) {
                DEBUG_VAR(edge->toString());
                DEBUG_VAR(edge->getVertexSrc()->toString());
                result = false;
            }
        } else {
            DEBUG_VAR(edge->toString());
            result = false;
        }
        if (edge->getVertexDst()) {
            if (edge->getVertexDst()->getEdgeIn() != edge) {
                DEBUG_VAR(edge->toString());
                DEBUG_VAR(edge->getVertexDst()->toString());
                result = false;
            }
        } else {
            DEBUG_VAR(edge->toString());
            result = false;
        }
    }
    return result;
}

double SphericalPolygon::getRadius() const {
    double result = 0.0;
    if (sphere_) {
#ifdef USE_CGAL
        result = CGAL::sqrt(sphere_->squared_radius());
#else
        result = sphere_->getRadius();
#endif
    }
    return result;
}

std::string SphericalPolygon::toString() const {
    std::stringstream sstr;
    sstr << "SphericalPolygon(";
    sstr << util::StringFactory::fromPointer(this) << ", ";
    Point3SPtr p_center = KernelFactory::createPoint3(sphere_);
    Vector3SPtr v_center = KernelFactory::createVector3(p_center);
    sstr << "<" << util::StringFactory::fromDouble((*v_center)[0]) << ", "
         << util::StringFactory::fromDouble((*v_center)[1]) << ", "
         << util::StringFactory::fromDouble((*v_center)[2]) << ">, ";
    sstr << "CircularVertices:" + util::StringFactory::fromInteger(vertices_.size()) + ", ";
    sstr << "CircularEdges:" + util::StringFactory::fromInteger(edges_.size()) + ",\n";
    std::list<CircularVertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        CircularVertexSPtr vertex = *it_v++;
        sstr << "\t" << vertex->toString() << "\n";
    }
    sstr << "\n";
    std::list<CircularEdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        CircularEdgeSPtr edge = *it_e++;
        sstr << "\t" << edge->toString() << "\n";
    }
    sstr << ")\n";
    return sstr.str();
}

} }
