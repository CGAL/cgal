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
 * @file   algo/3d/ConvexVertexSplitter.cpp
 * @author Gernot Walzl
 * @date   2013-02-01
 */

#include "algo/3d/ConvexVertexSplitter.h"

#include "debug.h"
#include "algo/3d/SelfIntersection.h"
#include "algo/3d/PolyhedronTransformation.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Polyhedron.h"
#include "db/3d/OBJFile.h"
#include "util/ptrs.h"
#include "util/Configuration.h"

#include <list>
#include <string>

namespace algo { namespace _3d {

ConvexVertexSplitter::ConvexVertexSplitter() {
    type_ = AbstractVertexSplitter::CONVEX_VERTEX_SPLITTER;
    optimization_ = -1;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        std::string s_optimization = config->getString(
                "algo_3d_ConvexVertexSplitter", "optimization");
        if (s_optimization.compare("max") == 0) {
            optimization_ = -1;
        } else if (s_optimization.compare("min") == 0) {
            optimization_ = 1;
        }
    }
}

ConvexVertexSplitter::~ConvexVertexSplitter() {
    // intentionally does nothing
}

ConvexVertexSplitterSPtr ConvexVertexSplitter::create() {
    ConvexVertexSplitterSPtr result =
            ConvexVertexSplitterSPtr(new ConvexVertexSplitter());
    return result;
}

int ConvexVertexSplitter::countConvexEdges(PolyhedronSPtr polyhedron) {
    int result = 0;
    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (!edge->isReflex()) {
            result++;
        }
    }
    return result;
}

PolyhedronSPtr ConvexVertexSplitter::splitVertex(VertexSPtr vertex) {
    std::cout << "\n> Splitting " << vertex->toString() << std::endl;

    PolyhedronSPtr polyhedron = vertex->getPolyhedron();
    if (vertex->degree() <= 3) {
        return polyhedron;
    }
    WriteLock l(polyhedron->mutex());
    vertex->sort();
    std::list<combi> combinations = generateAllCombinations(vertex->degree());
    std::cout << combinations.size() << " combinations" << std::endl;

    combi combi_opt;
    PolyhedronSPtr poly_opt;
    PolyhedronSPtr poly_opt_offset;
    int num_convex_edges = 0;
    int num_convex_edges_opt = 0;
    std::list<combi>::iterator it_combi = combinations.begin();
    while (it_combi != combinations.end()) {
        combi combination = *it_combi++;
        // std::cout << "-- Testing split-combination: " << combiToString(combination) << std::endl;

        // don't take it out of the loop
        PolyhedronSPtr poly_c = copyVertex(vertex);
        CGAL_assertion(bool(poly_c));
        CGAL_assertion(poly_c->facets().size() == vertex->facets().size());
        poly_c->initializeAllIDs();

        VertexSPtr vertex_c = poly_c->vertices().front();
        CombiVertexSplitter::splitVertex(vertex_c, combination);
        PolyhedronSPtr poly_c_offset = PolyhedronTransformation::shiftFacets(poly_c, -1.0);
        if (!poly_c_offset) {
            std::cerr << "Warning: failed to create offset of corner" << std::endl;
            continue;
        }

        // std::cout << "= Base Polyhedron\n" << poly_c->toString() << std::endl;
        // std::cout << "= Shifted Polyhedron\n" << poly_c_offset->toString() << std::endl;
        static int test_id = -1;
        ++test_id;
        // db::_3d::OBJFile::save("results/split_" + std::to_string(test_id) + ".obj", poly_c, false);
        // db::_3d::OBJFile::save("results/split_" + std::to_string(test_id) + "_offset.obj", poly_c_offset, false);
        // poly_c->dumpEdges("results/last_convex_split_base");
        // poly_c_offset->dumpEdges("results/last_convex_split_offset");

        if (!SelfIntersection::hasSelfIntersectingSurface(poly_c_offset)) {
            DEBUG_VAL("Valid split-combination found: " << combiToString(combination));
            num_convex_edges = countConvexEdges(poly_c_offset);
            DEBUG_VAR(num_convex_edges);
            if (!poly_opt) {
                combi_opt = combination;
                poly_opt = poly_c;
                poly_opt_offset = poly_c_offset;
                num_convex_edges_opt = num_convex_edges;
                continue;
            }
            if (optimization_ < 0) {
                if (num_convex_edges > num_convex_edges_opt) {
                    combi_opt = combination;
                    poly_opt = poly_c;
                    poly_opt_offset = poly_c_offset;
                    num_convex_edges_opt = num_convex_edges;
                }
            } else if (optimization_ > 0) {
                if (num_convex_edges < num_convex_edges_opt) {
                    combi_opt = combination;
                    poly_opt = poly_c;
                    poly_opt_offset = poly_c_offset;
                    num_convex_edges_opt = num_convex_edges;
                }
            }
        }
    }
    // CGAL_assertion(combi_opt != combi());
    DEBUG_VAL("Selected split-combination: " << combiToString(combi_opt));
    CombiVertexSplitter::apply(poly_opt, vertex);
    return polyhedron;
}

std::string ConvexVertexSplitter::toString() const {
    std::string result("ConvexVertexSplitter(");
    if (optimization_ == -1) {
        result += "max";
    } else if (optimization_ == 1) {
        result += "min";
    }
    result += ")";
    return result;
}

} }
