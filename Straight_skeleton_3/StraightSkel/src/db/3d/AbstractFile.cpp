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

#ifdef CGAL_SS3_DUMP_FILES
# include "db/3d/OBJFile.h"
#endif

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
        if(epsilon == 0.) {
            return (*(facet_l->plane()) == *(facet_r->plane())); // planes are normalized
        }

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

// @todo this doesn't really make sense if one or both of the facets is null
// (but we should never need it)
void AbstractFile::mergeFacets(EdgeSPtr edge,
                               FacetSPtr facet_into,
                               FacetSPtr facet_from,
                               PolyhedronSPtr polyhedron)
{
    CGAL_precondition(facet_into && facet_from);
    CGAL_precondition(facet_into != facet_from);
    CGAL_precondition(edge->getFacetL() == facet_from || edge->getFacetR() == facet_from);
    CGAL_precondition(edge->getFacetL() == facet_into || edge->getFacetR() == facet_into);

    CGAL_SS3_TRANSF_TRACE("Merging F" << facet_from->getID() << " into F" << facet_into->getID() <<
                          " Common edge E" << edge->getID() << " [V" << edge->getVertexSrc()->getID()
                                                            << " - V" << edge->getVertexDst()->getID() << "]");
    CGAL_SS3_TRANSF_TRACE("  FROM normal: " << *(KernelFactory::createVector3(facet_from->getPlane())));
    CGAL_SS3_TRANSF_TRACE("  INTO normal: " << *(KernelFactory::createVector3(facet_into->getPlane())));

    // remove the facet incidence info from the edge such that the facets are not deleted
    // when Polyhedron::removeEdge() is called
    facet_into->removeEdge(edge);
    facet_from->removeEdge(edge);

    polyhedron->removeEdge(edge);
    if (facet_into && facet_from && facet_into != facet_from) {
        facet_into->merge(facet_from);
        polyhedron->removeFacet(facet_from);
    }

    CGAL_postcondition(polyhedron->isConsistent());
}

void AbstractFile::mergeFacets(EdgeSPtr edge,
                               PolyhedronSPtr polyhedron)
{
    return mergeFacets(edge, edge->getFacetL(), edge->getFacetR(), polyhedron);
}

int AbstractFile::mergeCoplanarFacets(PolyhedronSPtr polyhedron,
                                      double epsilon)
{
    CGAL_SS3_TRANSF_TRACE("\nMerging coplanar faces with epsilon = " << epsilon);
    CGAL_SS3_TRANSF_TRACE("  initial face count: " << polyhedron->facets().size());

#ifdef CGAL_SS3_DUMP_FILES
    db::_3d::OBJFile::save("results/coplanar_merge_before.obj", polyhedron,
                           false /*do not triangulate*/);
#endif

    int result = 0;
    std::list<EdgeWPtr> edges_toremove;
    std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        EdgeSPtr edge = *it_e++;
        if (hasCoplanarFacets(edge, epsilon)) {
            edges_toremove.push_back(edge);
        }
    }

    CGAL_SS3_TRANSF_TRACE(edges_toremove.size() << " edges to merge");
    CGAL_SS3_TRANSF_TRACE_CODE(if (edges_toremove.size() > 0))
    CGAL_SS3_TRANSF_TRACE("Adjacent facets of the following edges are detected to be coplanar and will be merged.");

    std::list<EdgeWPtr>::iterator it_we = edges_toremove.begin();
    while (it_we != edges_toremove.end()) {
        EdgeSPtr edge = (it_we++)->lock();
        if (edge) {
            mergeFacets(edge, polyhedron);
            ++result;
        }
    }

    db::_3d::AbstractFile::sanitize(polyhedron);

    polyhedron->initializeAllIDs();

    CGAL_SS3_TRANSF_TRACE("  final face count: " << polyhedron->facets().size());

#ifdef CGAL_SS3_DUMP_FILES
    db::_3d::OBJFile::save("results/coplanar_merge_after.obj", polyhedron,
                           false /*do not triangulate*/);
#endif

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

int AbstractFile::removeFacetsDegLt3(PolyhedronSPtr polyhedron) {
    int result = 0;
    std::list<FacetSPtr> facets_tomerge;
    for (FacetSPtr facet : polyhedron->facets()) {
        if (facet->vertices().size() < 3) {
            facets_tomerge.push_back(facet);
        }
    }
    for (FacetSPtr facet : facets_tomerge) {
        // Facet could have grown from another merge, so check again
        if (facet->vertices().size() >= 3) {
            continue;
        }

        // Out of the two edges of the face, find the edge that is incident to the facet
        // that is the largest, and merge into that one
        EdgeSPtr best_edge = nullptr;
        FacetSPtr best_neighbor = nullptr;
        std::size_t best_size = 0;

        for (EdgeSPtr edge : facet->edges()) {
            FacetSPtr neighbor = nullptr;
            if (edge->getFacetL() == facet && edge->getFacetR() && edge->getFacetR() != facet) {
                neighbor = edge->getFacetR();
            } else if (edge->getFacetR() == facet && edge->getFacetL() && edge->getFacetL() != facet) {
                neighbor = edge->getFacetL();
            }
            if (neighbor && neighbor->vertices().size() > best_size) {
                best_edge = edge;
                best_neighbor = neighbor;
                best_size = neighbor->vertices().size();
            }
        }

        if (best_neighbor) {
            mergeFacets(best_edge, best_neighbor, facet, polyhedron);
        } else {
            // No neighbor, just remove
            for (EdgeSPtr edge : facet->edges()) {
                polyhedron->removeEdge(edge);
            }
        }
        ++result;
    }
    return result;
}

int AbstractFile::removeVerticesDegLt3(PolyhedronSPtr polyhedron) {
    CGAL_SS3_TRANSF_TRACE("Remove Vertices with degree < 3");
    CGAL_SS3_TRANSF_TRACE("  initial vertex count: " << polyhedron->vertices().size());

    int result = 0;
    std::list<VertexSPtr> vertices_toremove;
    std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
    while (it_v != polyhedron->vertices().end()) {
        VertexSPtr vertex = *it_v++;
        if (vertex->degree() < 3) {
            vertices_toremove.push_back(vertex);
            CGAL_SS3_TRANSF_TRACE("Enlist: V" << vertex->getID());
            for (FacetWPtr wf : vertex->facets()) {
                FacetSPtr facet = wf.lock();
                CGAL_SS3_TRANSF_TRACE("  Incident facet with: " << facet->vertices().size() << " vertices");
            }
        }
    }
    it_v = vertices_toremove.begin();
    while (it_v != vertices_toremove.end()) {
        VertexSPtr vertex = *it_v++;
        CGAL_SS3_TRANSF_TRACE("Removing " << vertex->toString());

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

            FacetSPtr fL = edge_src->getFacetL();
            FacetSPtr fR = edge_src->getFacetR();
            CGAL_assertion(fL != fR);

            // if there is any facet with degree 3, put it in fL (both could be with degree 3,
            // but then the swap does not change anything)
            //
            // It's "== 2" because we have already removed the vertex from its incident facets
            CGAL_assertion(!fL->hasVertex(vertex));
            CGAL_assertion(!fR->hasVertex(vertex));
            if (fR->vertices().size() == 2) {
                std::swap(fL, fR);
            }

            if (fL->vertices().size() == 2) {
                if (fR->vertices().size() == 2) {
                    CGAL_SS3_TRANSF_TRACE("Vertex is the apex of two facets of degree 3");
                    // both facets have degree 3, so remove everything (both facets, both edges,
                    // and one of the other edges + setting up incident facets properly)
                    VertexSPtr other_vertex = edge_src->other(vertex);
                    EdgeSPtr third_edge_1 = edge_src->next(other_vertex);
                    EdgeSPtr third_edge_2 = edge_src->prev(other_vertex);
                    CGAL_assertion(third_edge_1 != third_edge_2);
                    mergeFacets(third_edge_1, polyhedron);
                    mergeFacets(third_edge_2, polyhedron);

                    edge_dst->getFacetL()->removeEdge(edge_dst);
                    edge_dst->getFacetR()->removeEdge(edge_dst);
                    polyhedron->removeEdge(edge_dst);

                    if (edge_src->getVertexDst() == vertex) {
                        edge_src->replaceVertexDst(vertex_dst);
                    } else if (edge_src->getVertexSrc() == vertex) {
                        edge_src->replaceVertexSrc(vertex_dst);
                    }
                } else {
                    CGAL_SS3_TRANSF_TRACE("Vertex is the apex of one facet of degree 3");
                    // one facet has degree 3, so remove the vertex, the two incident edges, and the facet.
                    EdgeSPtr third_edge;
                    for (EdgeSPtr edge : fR->edges()) {
                        if (edge != edge_src && edge != edge_dst) {
                            third_edge = edge;
                            break;
                        }
                    }
                    CGAL_assertion(third_edge != nullptr);

                    fR->removeVertex(vertex);
                    fR->removeEdge(edge_src);
                    fR->removeEdge(edge_dst);
                    fL->removeEdge(third_edge);
                    if (third_edge->getFacetL() == fL) {
                        third_edge->setFacetL(fR);
                    } else {
                        third_edge->setFacetR(fR);
                    }
                }
            } else {
                edge_dst->getFacetL()->removeEdge(edge_dst);
                edge_dst->getFacetR()->removeEdge(edge_dst);
                polyhedron->removeEdge(edge_dst);

                if (edge_src->getVertexDst() == vertex) {
                    edge_src->replaceVertexDst(vertex_dst);
                } else if (edge_src->getVertexSrc() == vertex) {
                    edge_src->replaceVertexSrc(vertex_dst);
                }
            }
        }

        polyhedron->removeVertex(vertex);
        CGAL_postcondition(polyhedron->isConsistent());
        ++result;
    }

    CGAL_SS3_TRANSF_TRACE("  final vertex count: " << polyhedron->vertices().size());

    return result;
}

int AbstractFile::sanitize(PolyhedronSPtr polyhedron) {
    CGAL_SS3_TRANSF_TRACE("Sanitizing...");

    // - removeVerticesDegLt3 can create facets with fewer than 3 vertices
    // - removeFacetsDegLt3 removes facets with fewer than 3 vertices
    // so loop till nothing is done anymore
    int result = 0;
    for(;;) {
        std::size_t vlt3 = removeVerticesDegLt3(polyhedron);
        std::size_t flt3 = removeFacetsDegLt3(polyhedron);
        int partial = vlt3 + flt3;
        if (partial == 0) {
            break;
        }
        result += vlt3 + flt3;
    }
    return result;
}


} }
