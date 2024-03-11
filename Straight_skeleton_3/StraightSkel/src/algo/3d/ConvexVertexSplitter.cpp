/**
 * @file   algo/3d/ConvexVertexSplitter.cpp
 * @author Gernot Walzl
 * @date   2013-02-01
 */

#include "algo/3d/ConvexVertexSplitter.h"

#include "debug.h"
#include "algo/3d/SelfIntersection.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Polyhedron.h"
#include "util/ptrs.h"
#include "util/Configuration.h"
#include <list>

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
    PolyhedronSPtr polyhedron = vertex->getPolyhedron();
    if (vertex->degree() <= 3) {
        return polyhedron;
    }
    WriteLock l(polyhedron->mutex());
    vertex->sort();
    std::list<combi> combinations = generateAllCombinations(vertex->degree());
    combi combi_opt;
    PolyhedronSPtr poly_opt;
    PolyhedronSPtr poly_opt_offset;
    int num_convex_edges = 0;
    int num_convex_edges_opt = 0;
    std::list<combi>::iterator it_combi = combinations.begin();
    while (it_combi != combinations.end()) {
        combi combination = *it_combi++;
        PolyhedronSPtr poly_c = copyVertex(vertex);
        VertexSPtr vertex_c = poly_c->vertices().front();
        CombiVertexSplitter::splitVertex(vertex_c, combination);
        PolyhedronSPtr poly_c_offset = shiftFacets(poly_c, -1.0);
        if (!poly_c_offset) {
            continue;
        }
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
