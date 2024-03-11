/**
 * @file   algo/3d/VolumeVertexSplitter.cpp
 * @author Gernot Walzl
 * @date   2013-01-27
 */

#include "algo/3d/VolumeVertexSplitter.h"

#include "debug.h"
#include "algo/3d/SelfIntersection.h"
#include "algo/3d/KernelWrapper.h"
#include "data/2d/Vertex.h"
#include "data/2d/Edge.h"
#include "data/2d/Polygon.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "util/ptrs.h"
#include "util/Configuration.h"
#include <list>

namespace algo { namespace _3d {

VolumeVertexSplitter::VolumeVertexSplitter() {
    type_ = AbstractVertexSplitter::VOLUME_VERTEX_SPLITTER;
    optimization_ = 1;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        std::string s_optimization = config->getString(
                "algo_3d_VolumeVertexSplitter", "optimization");
        if (s_optimization.compare("max") == 0) {
            optimization_ = -1;
        } else if (s_optimization.compare("min") == 0) {
            optimization_ = 1;
        }
    }
}

VolumeVertexSplitter::~VolumeVertexSplitter() {
    // intentionally does nothing
}

VolumeVertexSplitterSPtr VolumeVertexSplitter::create() {
    VolumeVertexSplitterSPtr result =
            VolumeVertexSplitterSPtr(new VolumeVertexSplitter());
    return result;
}

double VolumeVertexSplitter::calcArea(FacetSPtr facet) {
    double result = 0.0;
    data::_2d::PolygonSPtr polygon_2d = facet->toPolygon();
    data::_2d::VertexSPtr vertex = polygon_2d->vertices().front();
    data::_2d::VertexSPtr vertex_begin;
    data::_2d::VertexSPtr vertex_prev;
    while (vertex != vertex_begin) {
        if (!vertex_begin) {
            vertex_begin = vertex;
        }
        if (vertex_prev && vertex) {
            result += vertex_prev->getX() * vertex->getY()
                    - vertex->getX() * vertex_prev->getY();
        }
        vertex_prev = vertex;
        vertex = vertex->next();
    }
    if (vertex_prev && vertex) {
        result += vertex_prev->getX() * vertex->getY()
                - vertex->getX() * vertex_prev->getY();
    }
    result *= 0.5;
    return result;
}

double VolumeVertexSplitter::calcSurfaceArea(PolyhedronSPtr polyhedron) {
    double result = 0.0;
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        result += calcArea(facet);
    }
    return result;
}


void VolumeVertexSplitter::closeFacets(PolyhedronSPtr polyhedron) {
    std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        FacetSPtr facet = *it_f++;
        VertexSPtr vertex_src = VertexSPtr();
        VertexSPtr vertex_dst = VertexSPtr();
        std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            if (edge->getVertexDst()->degree() < 3) {
                if (edge->getFacetL() == facet) {
                    vertex_src = edge->getVertexDst();
                } else if (edge->getFacetR() == facet) {
                    vertex_dst = edge->getVertexDst();
                }
            }
        }
        if (vertex_src && vertex_dst) {
            EdgeSPtr edge = Edge::create(vertex_src, vertex_dst);
            polyhedron->addEdge(edge);
            facet->addEdge(edge);
        } else {
            DEBUG_VAR(facet->toString());
        }
    }
}

int VolumeVertexSplitter::compareSurfaceAreas(PolyhedronSPtr polyhedron_1, PolyhedronSPtr polyhedron_2) {
    int result = 0;

    PolyhedronSPtr polyhedron_1_c = polyhedron_1->clone();
    PolyhedronSPtr polyhedron_2_c = polyhedron_2->clone();

    // set destinations of edges
    std::list<EdgeSPtr>::iterator it_e1 = polyhedron_1_c->edges().begin();
    while (it_e1 != polyhedron_1_c->edges().end()) {
        EdgeSPtr edge_1 = *it_e1++;
        VertexSPtr vertex_1_dst = edge_1->getVertexDst();
        if (vertex_1_dst->degree() > 1) {
            continue;
        }
        std::list<EdgeSPtr>::iterator it_e2 = polyhedron_2_c->edges().begin();
        while (it_e2 != polyhedron_2_c->edges().end()) {
            EdgeSPtr edge_2 = *it_e2++;
            VertexSPtr vertex_2_dst = edge_2->getVertexDst();
            if (vertex_2_dst->degree() > 1) {
                continue;
            }
            Vector3SPtr v_dir_1 = KernelFactory::createVector3(edge_1->line());
            Vector3SPtr v_dir_2 = KernelFactory::createVector3(edge_2->line());
            if (KernelWrapper::angle(v_dir_1, v_dir_2) < 0.0001) {
                // TODO: epsilon environment is not good
                // same edge
                Point3SPtr p_1_dst = vertex_1_dst->getPoint();
                Point3SPtr p_2_dst = vertex_2_dst->getPoint();
                Point3SPtr p_dst = p_1_dst;
                if (KernelWrapper::comparePoints(v_dir_1, p_1_dst, p_2_dst) < 0) {
                    p_dst = p_2_dst;
                }
                vertex_1_dst->setPoint(p_dst);
                vertex_2_dst->setPoint(p_dst);
                break;
            }
        }
    }

    closeFacets(polyhedron_1_c);
    closeFacets(polyhedron_2_c);
    double area_1 = calcSurfaceArea(polyhedron_1_c);
    double area_2 = calcSurfaceArea(polyhedron_2_c);
    if (area_1 < area_2) {
        result = -1;
    } else if (area_1 > area_2) {
        result = 1;
    }
    return result;
}


PolyhedronSPtr VolumeVertexSplitter::splitVertex(VertexSPtr vertex) {
    PolyhedronSPtr polyhedron = vertex->getPolyhedron();
    if (vertex->degree() <= 3) {
        return polyhedron;
    }
    WriteLock l(polyhedron->mutex());
    vertex->sort();
    std::list<combi> combinations = generateAllCombinations(vertex->degree());
    combi combi_opt_vol;
    PolyhedronSPtr poly_opt_vol;
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
            if (!poly_opt_vol) {
                combi_opt_vol = combination;
                poly_opt_vol = poly_c_offset;
                continue;
            }
            if (optimization_ > 0) {
                if (compareSurfaceAreas(poly_opt_vol, poly_c_offset) < 0) {
                    combi_opt_vol = combination;
                    poly_opt_vol = poly_c_offset;
                }
            } else if (optimization_ < 0) {
                if (compareSurfaceAreas(poly_opt_vol, poly_c_offset) > 0) {
                    combi_opt_vol = combination;
                    poly_opt_vol = poly_c_offset;
                }
            }
        }
    }
    DEBUG_VAL("Selected split-combination: " << combiToString(combi_opt_vol));
    CombiVertexSplitter::splitVertex(vertex, combi_opt_vol);
    return polyhedron;
}

std::string VolumeVertexSplitter::toString() const {
    std::string result("VolumeVertexSplitter(");
    if (optimization_ == -1) {
        result += "max";
    } else if (optimization_ == 1) {
        result += "min";
    }
    result += ")";
    return result;
}

} }
