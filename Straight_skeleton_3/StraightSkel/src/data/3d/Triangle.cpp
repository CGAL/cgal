/**
 * @file   data/3d/Triangle.cpp
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#include "data/3d/Triangle.h"

#include "data/3d/Vertex.h"
#include "data/3d/Facet.h"
#include "data/3d/KernelFactory.h"
#include "debug.h"
#include "util/StringFactory.h"

namespace data { namespace _3d {

Triangle::~Triangle() {
    // intentionally does nothing
}

TriangleSPtr Triangle::create(FacetSPtr facet, VertexSPtr vertices[3]) {
    TriangleSPtr result = TriangleSPtr(new Triangle());
    result->setFacet(facet);
    result->setVertices(vertices);
    facet->addTriangle(result);
    return result;
}

FacetSPtr Triangle::getFacet() const {
    DEBUG_WPTR(this->facet_);
    if (this->facet_.expired())
        return FacetSPtr();
    else
        return FacetSPtr(this->facet_);
}

void Triangle::setFacet(FacetSPtr facet) {
    this->facet_ = facet;
}

std::list<TriangleSPtr>::iterator Triangle::getFacetListIt() const {
    return this->facet_list_it_;
}

void Triangle::setFacetListIt(std::list<TriangleSPtr>::iterator list_it) {
    this->facet_list_it_ = list_it;
}

VertexSPtr Triangle::getVertex(unsigned int index) const {
    VertexSPtr result = VertexSPtr();
    if (0 <= index && index < 3) {
        result = vertices_[index];
    }
    DEBUG_SPTR(result);
    return result;
}

void Triangle::setVertices(VertexSPtr vertices[3]) {
    for (unsigned int i = 0; i < 3; i++) {
        this->vertices_[i] = vertices[i];
    }
}

Plane3SPtr Triangle::plane() const {
    return KernelFactory::createPlane3(vertices_[0]->getPoint(),
            vertices_[1]->getPoint(), vertices_[2]->getPoint());
}

int Triangle::getID() const {
    return this->id_;
}

void Triangle::setID(int id) {
    this->id_ = id;
}

std::string Triangle::toString() const {
    std::string result("Triangle(");
    for (unsigned int i = 0; i < 3; i++) {
        if (i > 0) {
            result += ", ";
        }
        if (vertices_[i]->getID() != -1) {
            result += util::StringFactory::fromInteger(vertices_[i]->getID());
        } else {
            result += util::StringFactory::fromPointer(vertices_[i].get());
        }
    }
    result += ")";
    return result;
}

} }
