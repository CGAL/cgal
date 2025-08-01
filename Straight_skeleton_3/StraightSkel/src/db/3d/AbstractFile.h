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
 * @file   db/3d/AbstractFile.h
 * @author Gernot Walzl
 * @date   2013-12-18
 */

#ifndef DB_3D_ABSTRACTFILE_H
#define DB_3D_ABSTRACTFILE_H

#include "debug.h"
#include "data/3d/ptrs.h"

#include "data/3d/KernelFactory.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

namespace db { namespace _3d {

using namespace data::_3d;

class AbstractFile {
public:
    virtual ~AbstractFile();

    static bool hasCoplanarFacets(EdgeSPtr edge, double epsilon);
    static void mergeFacets(EdgeSPtr edge, FacetSPtr facet_into, FacetSPtr facet_from, PolyhedronSPtr polyhedron);
    static void mergeFacets(EdgeSPtr edge, PolyhedronSPtr polyhedron);
    static int mergeCoplanarFacets(PolyhedronSPtr polyhedron, double epsilon);
    static int mergeCoplanarFacets(PolyhedronSPtr polyhedron);
    static int removeVerticesDegLt3(PolyhedronSPtr polyhedron);
    static int removeFacetsDegLt3(PolyhedronSPtr polyhedron);
    static int sanitize(PolyhedronSPtr polyhedron);

    // The tag is a template parameter because for debugging outputs,
    // we might want to tolerate intersections
    template <typename CDT2Tag>
    static auto constructFacetTriangulation(FacetSPtr facet)
    {
        using Itag = CDT2Tag;
        using PK = CGAL::Projection_traits_3<CGAL::K>;
        using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
        using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
        using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
        using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
        using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
        using PCDT_VH = typename PCDT::Vertex_handle;

        Vector3SPtr n = KernelFactory::createVector3(facet->plane());
        CGAL_precondition(*n != CGAL::NULL_VECTOR);

        PK projection_traits(*n);
        PCDT pcdt(projection_traits);

        std::map<VertexSPtr, PCDT_VH> face_vhs; // might have multiple vertices at the same position

        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            auto res = face_vhs.emplace(vertex, PCDT_VH());
            if(res.second) // first time seeing this point
            {
                PCDT_VH vh = pcdt.insert(*(vertex->getPoint()));
                res.first->second = vh;
                vh->info() = vertex;
            }
        }

        std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            VertexSPtr v0 = edge->src(facet);
            VertexSPtr v1 = edge->dst(facet);
            CGAL_assertion(*(v0->getPoint()) != *(v1->getPoint()));

            PCDT_VH vh0 = face_vhs.at(v0);
            PCDT_VH vh1 = face_vhs.at(v1);

            try {
                pcdt.insert_constraint(vh0, vh1);
            } catch(const typename PCDT::Intersection_of_constraints_exception&) {
                CGAL_SS3_TRANSF_TRACE("Error: Intersection of constraint w/ " << vh0->point() << " " << vh1->point());
                CGAL_SS3_TRANSF_TRACE(facet->toString());
                CGAL_assertion_msg(false, "Intersections in CDT2 are not allowed");
                return PCDT(projection_traits);
            }
        }

        return pcdt;
    }

protected:
    AbstractFile();
};

} }

#endif /* DB_3D_ABSTRACTFILE_H */
