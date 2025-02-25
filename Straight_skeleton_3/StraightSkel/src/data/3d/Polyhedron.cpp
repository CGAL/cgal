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
 * @file   data/3d/Polyhedron.cpp
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#include "data/3d/ptrs.h"
#include "data/3d/Polyhedron.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/skel/SkelFacetData.h"
#include "data/3d/Triangle.h"
#include "util/StringFactory.h"

#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <sstream>

namespace data { namespace _3d {

Polyhedron::Polyhedron() {
    this->id_ = -1;
}

Polyhedron::~Polyhedron() {
    facets_.clear();
    edges_.clear();
    vertices_.clear();
}

PolyhedronSPtr Polyhedron::create() {
    PolyhedronSPtr result = PolyhedronSPtr(new Polyhedron());
    return result;
}

PolyhedronSPtr Polyhedron::create(unsigned int num_facets, FacetSPtr facets[]) {
    PolyhedronSPtr result = PolyhedronSPtr(new Polyhedron());
    for (unsigned int i = 0; i < num_facets; i++) {
        result->addFacet(facets[i]);
    }
    return result;
}

PolyhedronSPtr Polyhedron::clone() const {
    std::map<VertexSPtr, VertexSPtr> vertices_c;
    std::map<EdgeSPtr, EdgeSPtr> edges_c;
    PolyhedronSPtr result = Polyhedron::create();
    result->setDescription(description_);
    std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        VertexSPtr vertex_c = vertex->clone();
        result->addVertex(vertex_c);
        vertices_c[vertex] = vertex_c;
    }
    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr src = vertices_c[edge->getVertexSrc()];
        VertexSPtr dst = vertices_c[edge->getVertexDst()];
        EdgeSPtr edge_c = Edge::create(src, dst);
        result->addEdge(edge_c);
        edges_c[edge] = edge_c;
    }
    std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        FacetSPtr facet = *it_f++;
        FacetSPtr facet_c = Facet::create();
        facet_c->setPlane(facet->getPlane());
        facet_c->setBasePlaneID(facet->getBasePlaneID());
        if (facet->hasData()) {
            skel::SkelFacetDataSPtr data = std::dynamic_pointer_cast<skel::SkelFacetData>(facet->getData());
            skel::SkelFacetDataSPtr data_c = skel::SkelFacetData::create(facet_c);
            data_c->setSpeed(data->getSpeed());
        }
        if (facet->getID() != -1) { // @todo remove this? Anyway the other IDs are not copied...
          facet_c->setID(facet->getID());
        }
        facet_c->cachedPlane_ = facet->cachedPlane_;
        facet_c->cachedSpeed_ = facet->cachedSpeed_;

        std::list<VertexSPtr>::const_iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            facet_c->addVertex(vertices_c[vertex]);
        }
        std::list<EdgeSPtr>::const_iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            EdgeSPtr edge_c = edges_c[edge];
            if (edge->getFacetL() == facet) {
                edge_c->setFacetL(facet_c);
            }
            if (edge->getFacetR() == facet) {
                edge_c->setFacetR(facet_c);
            }
            facet_c->addEdge(edge_c);
        }
        std::list<TriangleSPtr>::const_iterator it_t = facet->triangles().begin();
        while (it_t != facet->triangles().end()) {
            TriangleSPtr triangle = *it_t++;
            VertexSPtr vertices_t_c[3];
            for (unsigned int i = 0; i < 3; i++) {
                vertices_t_c[i] = vertices_c[triangle->getVertex(i)];
            }
            Triangle::create(facet_c, vertices_t_c);
        }
        result->addFacet(facet_c);
    }
    return result;
}

void Polyhedron::addVertex(VertexSPtr vertex) {
    std::list<VertexSPtr>::iterator it = vertices_.insert(vertices_.end(), vertex);
    vertex->setPolyhedron(shared_from_this());
    vertex->setPolyhedronListIt(it);
}

bool Polyhedron::removeVertex(VertexSPtr vertex) {
    bool result = false;
    if (vertex->getPolyhedron() == shared_from_this()) {
        vertices_.erase(vertex->getPolyhedronListIt());
        vertex->setPolyhedronListIt(std::list<VertexSPtr>::iterator());
        vertex->setPolyhedron(PolyhedronSPtr());
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet = FacetSPtr(facet_wptr);
                this->removeFacet(facet);
            }
        }
        std::list<EdgeWPtr>::iterator it_e = vertex->edges().begin();
        while (it_e != vertex->edges().end()) {
            EdgeWPtr edge_wptr = *it_e++;
            if (!edge_wptr.expired()) {
                EdgeSPtr edge = EdgeSPtr(edge_wptr);
                this->removeEdge(edge);
            }
        }
        result = true;
    }
    return result;
}

VertexSPtr Polyhedron::findVertex(VertexSPtr needle) {
    VertexSPtr result = VertexSPtr();
    std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        if (vertex->getPoint() == needle->getPoint()) {
            result = vertex;
            break;
        }
    }
    return result;
}

void Polyhedron::addEdge(EdgeSPtr edge) {
    std::list<EdgeSPtr>::iterator it = edges_.insert(edges_.end(), edge);
    edge->setPolyhedron(shared_from_this());
    edge->setPolyhedronListIt(it);
    VertexSPtr vertex = edge->getVertexSrc();
    if (vertex->getPolyhedron() != shared_from_this()) {
        this->addVertex(vertex);
    }
    vertex = edge->getVertexDst();
    if (vertex->getPolyhedron() != shared_from_this()) {
        this->addVertex(vertex);
    }
}

bool Polyhedron::removeEdge(EdgeSPtr edge) {
    bool result = false;
    if (edge->getPolyhedron() == shared_from_this()) {
        edges_.erase(edge->getPolyhedronListIt());
        edge->setPolyhedronListIt(std::list<EdgeSPtr>::iterator());
        edge->setPolyhedron(PolyhedronSPtr());
        FacetSPtr facet = edge->getFacetL();
        if (facet) {
            this->removeFacet(facet);
        }
        facet = edge->getFacetR();
        if (facet) {
            this->removeFacet(facet);
        }
        edge->getVertexSrc()->removeEdge(edge);
        edge->getVertexDst()->removeEdge(edge);
        result = true;
    }
    return result;
}

EdgeSPtr Polyhedron::findEdge(EdgeSPtr needle) {
    EdgeSPtr result = EdgeSPtr();
    std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getVertexSrc()->getPoint() ==
                needle->getVertexSrc()->getPoint() &&
                edge->getVertexDst()->getPoint() ==
                needle->getVertexDst()->getPoint()) {
            result = edge;
            break;
        }
        if (edge->getVertexSrc()->getPoint() ==
                needle->getVertexDst()->getPoint() &&
                edge->getVertexDst()->getPoint() ==
                needle->getVertexSrc()->getPoint()) {
            result = edge;
            break;
        }
    }
    return result;
}

void Polyhedron::addFacet(FacetSPtr facet) {
    std::list<FacetSPtr>::iterator it_f = facets_.insert(facets_.end(), facet);
    facet->setPolyhedronListIt(it_f);
    facet->setPolyhedron(shared_from_this());
    // add content of facet
    std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
    while (it_v != facet->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (vertex->getPolyhedron() != shared_from_this()) {
            this->addVertex(vertex);
        }
    }
    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getPolyhedron() != shared_from_this()) {
            this->addEdge(edge);
        }
    }
}

bool Polyhedron::removeFacet(FacetSPtr facet) {
    bool result = false;
    if (facet->getPolyhedron() == shared_from_this()) {
        facets_.erase(facet->getPolyhedronListIt());
        facet->setPolyhedron(PolyhedronSPtr());
        facet->setPolyhedronListIt(std::list<FacetSPtr>::iterator());
        std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            facet->removeEdge(edge);  // clears list iterators
        }
        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            facet->removeVertex(vertex);
        }
        result = true;
    }
    return result;
}

void Polyhedron::initPlanes() {
    std::list<FacetSPtr>::iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        FacetSPtr facet = *it_f++;
        facet->initPlane();
    }
}

void Polyhedron::clearData() {
    std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        vertex->setData(VertexDataSPtr());
    }
    std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        edge->setData(EdgeDataSPtr());
    }
    std::list<FacetSPtr>::iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        FacetSPtr facet = *it_f++;
        facet->setData(FacetDataSPtr());
    }
}

SharedMutex& Polyhedron::mutex() {
    return this->mutex_;
}

std::list<VertexSPtr>& Polyhedron::vertices() {
    return this->vertices_;
}

std::list<EdgeSPtr>& Polyhedron::edges() {
    return this->edges_;
}

std::list<FacetSPtr>& Polyhedron::facets() {
    return this->facets_;
}

VertexSPtr Polyhedron::ith_vertex(const std::size_t i)
{
    CGAL_assertion(i < this->vertices_.size());
    return *(std::next(std::cbegin(this->vertices_), i));
}
EdgeSPtr Polyhedron::ith_edge(const std::size_t i)
{
    CGAL_assertion(i < this->edges_.size());
    return *(std::next(std::cbegin(this->edges_), i));
}
FacetSPtr Polyhedron::ith_facet(const std::size_t i)
{
    CGAL_assertion(i < this->facets_.size());
    return *(std::next(std::cbegin(this->facets_), i));
}

bool Polyhedron::isConsistent() const {
    bool result = true;

    std::list<VertexSPtr>::const_iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        if (vertex->getPolyhedron() != shared_from_this()) {
            DEBUG_VAR(vertex->toString());
            result = false;
            break;
        }
        std::list<EdgeWPtr>::const_iterator it_e = vertex->edges().begin();
        unsigned int vne = 0;
        while (it_e != vertex->edges().end()) {
            EdgeWPtr edge_wptr = *it_e++;
            if (edge_wptr.expired()) {
                DEBUG_VAR(vertex->toString());
            } else {
                EdgeSPtr edge = EdgeSPtr(edge_wptr);
                if (vertex != edge->getVertexSrc() &&
                        vertex != edge->getVertexDst()) {
                    DEBUG_VAR(vertex->toString());
                    DEBUG_VAR(edge->toString());
                    result = false;
                    break;
                }

                if (edge->getVertexSrc() == edge->getVertexDst())
                {
                    DEBUG_VAR(edge->getVertexSrc());
                    result = false;
                    break;
                }

                // if (edge->getVertexSrc()->getPoint() == edge->getVertexDst()->getPoint())
                // {
                //     std::cerr << "- Degenerate edge!" << std::endl;
                //     std::cerr << edge->getVertexSrc()->toString() << std::endl;
                //     std::cerr << edge->getVertexDst()->toString() << std::endl;
                // }

                ++vne;
            }
        }

        // if (vne != 3) {
        //     std::cerr << "Warning: Vertex does not have 3 incident edges: " << vne << std::endl;
        //     std::cerr << vertex->toString() << std::endl;
        // }

        std::list<FacetWPtr>::const_iterator it_f = vertex->facets().begin();
        unsigned int vnf = 0;
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (facet_wptr.expired()) {
                DEBUG_VAR(vertex->toString());
            } else {
                FacetSPtr facet = FacetSPtr(facet_wptr);
                if (!facet->containsVertex(vertex)) {
                    DEBUG_VAR(vertex->toString());
                    DEBUG_VAR(facet->toString());
                    result = false;
                    break;
                }
                ++vnf;
            }
        }

        // if (vnf != 3) {
        //     std::cerr << "Warning: Vertex with not 3 incident faces? " << vnf << std::endl;
        //     std::cerr << vertex->toString() << std::endl;
        // }
    }

    std::list<EdgeSPtr>::const_iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        if (edge->getPolyhedron() != shared_from_this()) {
            DEBUG_VAR(edge->toString());
            result = false;
            break;
        }
        if (!edge->getVertexSrc()->containsEdge(edge)) {
            DEBUG_VAR(edge->toString());
            DEBUG_VAR(edge->getVertexSrc()->toString());
            result = false;
            break;
        }
        if (!edge->getVertexDst()->containsEdge(edge)) {
            DEBUG_VAR(edge->toString());
            DEBUG_VAR(edge->getVertexDst()->toString());
            result = false;
            break;
        }
        EdgeWPtr edge_wptr;
        edge_wptr = *(edge->getVertexSrcListIt());
        if (EdgeSPtr(edge_wptr) != edge) {
            DEBUG_VAR(edge->toString());
        }
        edge_wptr = *(edge->getVertexDstListIt());
        if (EdgeSPtr(edge_wptr) != edge) {
            DEBUG_VAR(edge->toString());
        }
        if (edge->getFacetL()) {
            if (!edge->getFacetL()->containsEdge(edge)) {
                DEBUG_VAR(edge->toString());
                DEBUG_VAR(edge->getFacetL()->toString());
                result = false;
                break;
            }
            edge_wptr = *(edge->getFacetLListIt());
            if (EdgeSPtr(edge_wptr) != edge) {
                DEBUG_VAR(edge->toString());
            }
        } else {
            std::cerr << "no left facet?" << std::endl;
            DEBUG_VAR(edge->toString());
            result = false;
        }
        if (edge->getFacetR()) {
            if (!edge->getFacetR()->containsEdge(edge)) {
                DEBUG_VAR(edge->toString());
                DEBUG_VAR(edge->getFacetR()->toString());
                result = false;
                break;
            }
            edge_wptr = *(edge->getFacetRListIt());
            if (EdgeSPtr(edge_wptr) != edge) {
                DEBUG_VAR(edge->toString());
            }
        } else {
            std::cerr << "no right facet?" << std::endl;
            DEBUG_VAR(edge->toString());
            result = false;
        }
    }

    std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        FacetSPtr facet = *it_f++;
        if (facet->getPolyhedron() != shared_from_this()) {
            DEBUG_VAR(facet->toString());
            result = false;
            break;
        }
        std::list<VertexSPtr>::const_iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            if (!vertex->containsFacet(facet)) {
                DEBUG_VAR(facet->toString());
                DEBUG_VAR(vertex->toString());
                result = false;
                break;
            }
        }
        std::list<EdgeSPtr>::const_iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            if (edge->getFacetL() != facet && edge->getFacetR() != facet) {
                DEBUG_VAR(facet->toString());
                DEBUG_VAR(edge->toString());
                result = false;
                break;
            }
            if (!facet->containsVertex(edge->getVertexSrc()) ||
                    !facet->containsVertex(edge->getVertexDst())) {
                DEBUG_VAR(facet->toString());
                DEBUG_VAR(edge->toString());
                result = false;
                break;
            }
        }
    }

    return result;
}

void Polyhedron::clear() {
    std::list<FacetSPtr>::iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        FacetSPtr facet = *it_f++;
        this->removeFacet(facet);
    }
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

bool Polyhedron::empty() {
  return vertices_.empty() && edges_.empty() && facets_.empty();
}

int Polyhedron::getID() const {
    return this->id_;
}

void Polyhedron::setID(int id) {
    this->id_ = id;
}

void Polyhedron::initializeAllIDs() {
    resetAllIDs();

    unsigned int vertex_id = 0;
    unsigned int edge_id = 0;
    unsigned int facet_id = 0;

    std::list<FacetSPtr>::iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        FacetSPtr facet = *it_f++;
        std::list<TriangleSPtr>::iterator it_t = facet->triangles().begin();
        while (it_t != facet->triangles().end()) {
            TriangleSPtr triangle = *it_t++;
            triangle->setID(-1);
        }
        facet->setID(facet_id++);
    }
    std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        edge->setID(edge_id++);
    }
    std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        vertex->setID(vertex_id++);
    }
    setID(-1);
}

void Polyhedron::resetAllIDs() {
    std::list<FacetSPtr>::iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        FacetSPtr facet = *it_f++;
        std::list<TriangleSPtr>::iterator it_t = facet->triangles().begin();
        while (it_t != facet->triangles().end()) {
            TriangleSPtr triangle = *it_t++;
            triangle->setID(-1);
        }
        facet->setID(-1);
    }
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

std::string Polyhedron::getDescription() const {
    return this->description_;
}

void Polyhedron::setDescription(const std::string& desc) {
    this->description_ = desc;
}

void Polyhedron::appendDescription(const std::string& desc) {
    this->description_.append(desc);
}

std::string Polyhedron::toString() const {
    std::stringstream sstr;
    sstr << "Polyhedron(";
    if (id_ != -1) {
        sstr << "id=" << util::StringFactory::fromInteger(id_) << ", ";
    } else {
        // sstr << util::StringFactory::fromPointer(this) << ", ";
    }
    sstr << "Vertices:" + util::StringFactory::fromInteger(vertices_.size()) + ", ";
    sstr << "Edges:" + util::StringFactory::fromInteger(edges_.size()) + ", ";
    sstr << "Facets:" + util::StringFactory::fromInteger(facets_.size()) + ",\n";

    std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        FacetSPtr facet = *it_f++;
        sstr << facet->toString() << "\n";
    }
    sstr << ")";

  /*std::list<EdgeSPtr>::iterator it_e = edges_.begin();
    while (it_e != edges_.end()) {
        EdgeSPtr edge = *it_e++;
        sstr << edge->toString() << "\n";
    }
    std::list<VertexSPtr>::iterator it_v = vertices_.begin();
    while (it_v != vertices_.end()) {
        VertexSPtr vertex = *it_v++;
        sstr << vertex->toString() << "\n";
    }*/
    return sstr.str();
}

void Polyhedron::dumpEdges(const std::string& base_filename) const
{
    unsigned int fi = 0;
    std::list<FacetSPtr>::const_iterator it_f = facets_.begin();
    while (it_f != facets_.end()) {
        FacetSPtr facet = *it_f++;

        std::ofstream edge_out(base_filename + "_" + std::to_string(fi++) + ".polylines.txt");
        edge_out.precision(17);

        std::list<EdgeSPtr>::const_iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            edge_out << "2 " << *(edge->getVertexSrc()->getPoint()) << " "
                             << *(edge->getVertexDst()->getPoint()) << std::endl;
        }
    }
}

} }
