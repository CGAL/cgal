/**
 * @file   data/2d/Polygon.cpp
 * @author Gernot Walzl
 * @date   2011-11-22
 */

#include "data/2d/Polygon.h"

#include "data/2d/Edge.h"
#include "data/2d/Vertex.h"
#include "debug.h"
#include "util/StringFactory.h"
#include <sstream>

namespace data { namespace _2d {

Polygon::Polygon() {
    this->id_ = -1;
}

Polygon::~Polygon() {
    edges_.clear();
    vertices_.clear();
}

PolygonSPtr Polygon::create() {
    PolygonSPtr result = PolygonSPtr(new Polygon());
    return result;
}

PolygonSPtr Polygon::create(unsigned int num_edges, EdgeSPtr edges[]) {
    PolygonSPtr result = PolygonSPtr(new Polygon());
    for (unsigned int i = 0; i < num_edges; i++) {
        result->addEdge(edges[i]);
    }
    result->sortEdges();
    return result;
}

void Polygon::addVertex(VertexSPtr vertex) {
    std::list<VertexSPtr>::iterator it = vertices_.insert(vertices_.end(), vertex);
    vertex->setPolygon(shared_from_this());
    vertex->setListIt(it);
}

bool Polygon::removeVertex(VertexSPtr vertex) {
    bool result = false;
    if (vertex->getPolygon() == shared_from_this()) {
        vertices_.erase(vertex->getListIt());
        vertex->setPolygon(PolygonSPtr());
        vertex->setListIt(std::list<VertexSPtr>::iterator());
        EdgeSPtr edge = vertex->getEdgeIn();
        if (edge) {
            this->removeEdge(edge);
            vertex->setEdgeIn(EdgeSPtr());
        }
        edge = vertex->getEdgeOut();
        if (edge) {
            this->removeEdge(edge);
            vertex->setEdgeOut(EdgeSPtr());
        }
        result = true;
    }
    return result;
}

void Polygon::addEdge(EdgeSPtr edge) {
    std::list<EdgeSPtr>::iterator it = edges_.insert(edges_.end(), edge);
    edge->setPolygon(shared_from_this());
    edge->setListIt(it);
    VertexSPtr vertex = edge->getVertexSrc();
    if (vertex->getPolygon() != shared_from_this()) {
        this->addVertex(vertex);
    }
    vertex = edge->getVertexDst();
    if (vertex->getPolygon() != shared_from_this()) {
        this->addVertex(vertex);
    }
}

bool Polygon::removeEdge(EdgeSPtr edge) {
    bool result = false;
    if (edge->getPolygon() == shared_from_this()) {
        edges_.erase(edge->getListIt());
        edge->setPolygon(PolygonSPtr());
        edge->setListIt(std::list<EdgeSPtr>::iterator());
        VertexSPtr vertex = edge->getVertexSrc();
        if (vertex->getEdgeOut() == edge) {
            vertex->setEdgeOut(EdgeSPtr());
        }
        vertex = edge->getVertexDst();
        if (vertex->getEdgeIn() == edge) {
            vertex->setEdgeIn(EdgeSPtr());
        }
        result = true;
    }
    return result;
}

void Polygon::sortEdges() {
    std::list<EdgeSPtr> temp;
    while (edges_.size() > 0) {
        EdgeSPtr first = EdgeSPtr();
        EdgeSPtr edge = edges_.front();
        while (edge != first) {
            if (!first) {
                first = edge;
            }
            edges_.erase(edge->getListIt());
            temp.push_back(edge);
            edge = edge->getVertexDst()->getEdgeOut();
        }
    }
    edges_.clear();
    std::list<EdgeSPtr>::iterator it_e = temp.begin();
    while (it_e != temp.end()) {
        EdgeSPtr edge = *it_e++;
        edge->setListIt(edges_.insert(edges_.end(), edge));
    }
}

void Polygon::clear() {
    std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        this->removeEdge(edge);
    }
    std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        this->removeVertex(vertex);
    }
}

SharedMutex& Polygon::mutex() {
    return this->mutex_;
}

std::list<VertexSPtr>& Polygon::vertices() {
    return this->vertices_;
}

std::list<EdgeSPtr>& Polygon::edges() {
    return this->edges_;
}

int Polygon::getID() const {
    return this->id_;
}

void Polygon::setID(int id) {
    this->id_ = id;
}

void Polygon::resetAllIDs() {
    std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        edge->setID(-1);
    }
    std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        vertex->setID(-1);
    }
    setID(-1);
}

bool Polygon::isConsistent() const {
    bool result = true;
    std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
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
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
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

int Polygon::countReflex() const {
    int result = 0;
    std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        if (vertex->isReflex()) {
            result += 1;
        }
    }
    return result;
}

int Polygon::countHoles() const {
    int result = -1;
    EdgeSPtr edge;
    EdgeSPtr edge_prev;
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        edge = *it_e++;
        if (edge->getVertexSrc()->getEdgeIn() != edge_prev) {
            result += 1;
        }
        edge_prev = edge;
    }
    return result;
}

std::string Polygon::getDescription() const {
    return this->description_;
}

void Polygon::setDescription(const std::string& desc) {
    this->description_ = desc;
}

void Polygon::appendDescription(const std::string& desc) {
    this->description_.append(desc);
}

std::string Polygon::toString() const {
    std::stringstream sstr;
    sstr << "Polygon(";
    if (id_ != -1) {
        sstr << "id=" << util::StringFactory::fromInteger(id_) << ", ";
    } else {
        sstr << util::StringFactory::fromPointer(this) << ", ";
    }
    sstr << "Vertices:" + util::StringFactory::fromInteger(vertices_.size()) + ", ";
    sstr << "Edges:" + util::StringFactory::fromInteger(edges_.size()) + ",\n";
    std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        sstr << "\t" << vertex->toString() << "\n";
    }
    sstr << "\n";
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
        while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        sstr << "\t" << edge->toString() << "\n";
    }
    sstr << ")\n";
    return sstr.str();
}

} }
