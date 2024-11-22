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
 * @file   db/3d/AbstractFile.cpp
 * @author Gernot Walzl
 * @date   2013-12-18
 */

#include "db/3d/AbstractFile.h"

#include "util/Configuration.h"
#include "debug.h"

#include "data/3d/KernelFactory.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "db/3d/OBJFile.h"

#include <cmath>
#include <list>

namespace db { namespace _3d {

AbstractFile::AbstractFile() {
    // intentionally does nothing
}

AbstractFile::~AbstractFile() {
    // intentionally does nothing
}

bool AbstractFile::hasCoplanarFacets(EdgeSPtr edge, double epsilon) {
    bool result = false;
    FacetSPtr facet_l = edge->getFacetL();
    FacetSPtr facet_r = edge->getFacetR();
    if (facet_l && facet_r) {
        if(epsilon == 0)
          return (*(facet_l->plane()) == *(facet_r->plane())); // planes are normalized

        Vector3SPtr normal_l = KernelFactory::createVector3(facet_l->plane());
        Vector3SPtr normal_r = KernelFactory::createVector3(facet_r->plane());
        CGAL::FT length_l = 0.0;
        CGAL::FT length_r = 0.0;
        for (unsigned int i = 0; i < 3; i++) {
            length_l += (*normal_l)[i] * (*normal_l)[i];
            length_r += (*normal_r)[i] * (*normal_r)[i];
        }
        // @todo tolerate this sqrt for now because it does not matter for robustness
        length_l = CGAL::sqrt_with_warning(length_l);
        length_r = CGAL::sqrt_with_warning(length_r);
        CGAL::FT diff = 0.0;
        CGAL::FT diff_sq_length = 0.0;
        for (unsigned int i = 0; i < 3; i++) {
            diff = ((*normal_l)[i]/length_l) - ((*normal_r)[i]/length_r);
            diff_sq_length += diff*diff;
        }
        result = (diff_sq_length < CGAL::square(epsilon));
    }
    return result;
}

int AbstractFile::mergeCoplanarFacets(PolyhedronSPtr polyhedron,
                                      double epsilon)
{
    std::cout << "\n> Merging coplanar faces with epsilon = " << epsilon << std::endl;
    std::cout << "  initial face count: " << polyhedron->facets().size() << std::endl;

    db::_3d::OBJFile::save("results/coplanar_merge_before.obj", polyhedron, false /*do not triangulate*/);

    int result = 0;
    std::list<EdgeSPtr> edges_toremove;
    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (hasCoplanarFacets(edge, epsilon)) {
            edges_toremove.push_back(edge);
        }
    }
    std::cout << edges_toremove.size() << " edges to merge" << std::endl;
    if (edges_toremove.size() > 0) {
        DEBUG_PRINT("Adjacent facets of the following edges are detected to be coplanar and will be merged.");
    }
    it_e = edges_toremove.begin();
    while (it_e != edges_toremove.end()) {
        EdgeSPtr edge = *it_e++;
        DEBUG_VAR(edge->toString());
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        if (facet_l) {
            facet_l->removeEdge(edge);
        }
        if (facet_r) {
            facet_r->removeEdge(edge);
        }
        polyhedron->removeEdge(edge);
        if (facet_l && facet_r &&
                facet_l != facet_r) {
            facet_l->merge(facet_r);
            polyhedron->removeFacet(facet_r);
        }
        result++;
    }
    std::cout << "  final face count: " << polyhedron->facets().size() << std::endl;
    db::_3d::OBJFile::save("results/coplanar_merge_after.obj", polyhedron, false /*do not triangulate*/);
    return result;
}

int AbstractFile::removeVerticesDegLt3(PolyhedronSPtr polyhedron) {
    int result = 0;
    std::list<VertexSPtr> vertices_toremove;
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        std::cout << "Considering " << vertex->toString() << std::endl;
        if (vertex->degree() < 3) {
            vertices_toremove.push_back(vertex);
        }
    }
    it_v = vertices_toremove.begin();
    while (it_v != vertices_toremove.end()) {
        VertexSPtr vertex = *it_v++;
        std::cout << "Removing " << vertex->toString() << std::endl;
        std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
        while (it_f != vertex->facets().end()) {
            FacetWPtr facet_wptr = *it_f++;
            if (!facet_wptr.expired()) {
                FacetSPtr facet(facet_wptr);
                facet->removeVertex(vertex);
            }
        }

        if (vertex->degree() == 2) {
            EdgeSPtr edge_src = vertex->firstEdge();
            EdgeSPtr edge_dst = edge_src->next(vertex);
            VertexSPtr vertex_src = edge_src->getVertexSrc();
            if (vertex_src == vertex) {
                vertex_src = edge_src->getVertexDst();
            }
            VertexSPtr vertex_dst = edge_dst->getVertexSrc();
            if (vertex_dst == vertex) {
                vertex_dst = edge_dst->getVertexDst();
            }
            edge_dst->getFacetL()->removeEdge(edge_dst);
            edge_dst->getFacetR()->removeEdge(edge_dst);
            polyhedron->removeEdge(edge_dst);
            if (edge_src->getVertexDst() == vertex) {
                edge_src->replaceVertexDst(vertex_dst);
            } else if (edge_src->getVertexSrc() == vertex) {
                edge_src->replaceVertexSrc(vertex_dst);
            }
        }
        polyhedron->removeVertex(vertex);
        result++;
    }
    return result;
}

int AbstractFile::mergeCoplanarFacets(PolyhedronSPtr polyhedron)
{
    double epsilon = 0.0;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        std::string section("main");
        std::string key("epsilon_coplanarity");
        if (config->contains(section, key)) {
            epsilon = config->getDouble(section, key);
        }
    }

    return mergeCoplanarFacets(polyhedron, epsilon);
}

} }
