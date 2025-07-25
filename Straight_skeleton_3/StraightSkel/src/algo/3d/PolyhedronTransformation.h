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
 * @file   algo/3d/PolyhedronTransformation.h
 * @author Gernot Walzl
 * @date   2012-09-01
 */

#ifndef ALGO_3D_POLYHEDRONTRANSFORMATION_H
#define ALGO_3D_POLYHEDRONTRANSFORMATION_H

#include "debug.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "db/3d/Surface_meshIO.h"
#include "util/Configuration.h"

#include <CGAL/Surface_mesh.h>

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

class PolyhedronTransformation {
public:
    virtual ~PolyhedronTransformation();

    static Point3SPtr boundingBoxMin(PolyhedronSPtr polyhedron);
    static Point3SPtr boundingBoxMax(PolyhedronSPtr polyhedron);

    static bool isInsideBox(PolyhedronSPtr polyhedron,
                            Point3SPtr p_box_min, Point3SPtr p_box_max);

    static void translate(PolyhedronSPtr polyhedron, Vector3SPtr t);
    static void scale(PolyhedronSPtr polyhedron, Vector3SPtr s);

    static void translateNscale(PolyhedronSPtr polyhedron,
                                Point3SPtr p_box_min, Point3SPtr p_box_max);

    static void truncatePrecision(PolyhedronSPtr polyhedron);

    /**
     * updates the position of the vertex of a polyhedron, computed from the planes of
     * its incident faces.
     */
    static bool resetPoint(VertexSPtr vertex, const std::array<Plane3SPtr, 3>& planes);
    static bool resetPoint(VertexSPtr vertex);

    /**
     * updates the positions of the vertices of a polyhedron, computed from the planes of
     * their incident faces.
     */
    static bool resetPoints(PolyhedronSPtr polyhedron);

    /**
     * returns the shifted position of the vertex of a polyhedron
     */
    static Point3SPtr shiftPoint(VertexSPtr vertex, const CGAL::FT& offset);

    /**
     * returns the shifted position of the edge of a polyhedron
     */
    static Segment3SPtr shiftEdge(EdgeSPtr edge, const CGAL::FT& offset);

    /**
     * returns the shifted position of the facet of a polyhedron
     */
    static Plane3SPtr shiftPlane(FacetSPtr vertex, const CGAL::FT& offset);

    /**
     * Offsets the polyhedron `polyhedron`
     * Negative offset points to the interior of the polyhedron.
     */
    static void shiftFacetsInPlace(PolyhedronSPtr polyhedron,
                                   const CGAL::FT& offset,
                                   const bool recompute_positions = true);

    /**
     * Creates an offset polyhedron.
     * Negative offset points to the interior of the polyhedron.
     */
    static PolyhedronSPtr shiftFacets(PolyhedronSPtr polyhedron,
                                      const CGAL::FT& offset,
                                      const bool recompute_positions = true);


    /**
     * Triangulate the facet 'f' and returns all newly created triangle facets
     */
    static std::pair<std::list<VertexSPtr>, std::list<FacetSPtr> > triangulate(FacetSPtr f,
                                            PolyhedronSPtr polyhedron);

    /**
     * Normalize facet planes
    */
    static void normalizeFacetPlanes(PolyhedronSPtr polyhedron);

    /**
     * Checking for degeneracies
     */
    static bool doAll2PlanesIntersect(PolyhedronSPtr polyhedron);
    static bool doAll3PlanesIntersect(PolyhedronSPtr polyhedron);

    static bool areAllVerticesDegree3(PolyhedronSPtr polyhedron);

    /**
     * Check that all faces have at most two high degree vertices
     */
    static bool isTiltCompatible(PolyhedronSPtr polyhedron);

    static void randMovePoints(PolyhedronSPtr polyhedron);
    static void randTiltPlanes(PolyhedronSPtr polyhedron);
    static void randTiltPlanesv3(PolyhedronSPtr polyhedron);

    template <typename NamedParameters = CGAL::parameters::Default_named_parameters>
    static PolyhedronSPtr convert(const CGAL::Surface_mesh<Point3>& sm,
                                  const NamedParameters& np = CGAL::parameters::default_values());

protected:
    static std::array<double, 3> randVec(double min, double max);
    PolyhedronTransformation();
};

template <typename NamedParameters>
PolyhedronSPtr
PolyhedronTransformation::
convert(const CGAL::Surface_mesh<Point3>& sm,
        const NamedParameters& np)
{
    CGAL_SS3_TRANSF_TRACE("Converting mesh...");

    bool merge_faces = false;

    util::ConfigurationSPtr config = util::Configuration::getInstance();
    std::string section("main");
    if (config->isLoaded() &&
        config->contains(section, "merge_coplanar_faces") &&
        config->getBool(section, "merge_coplanar_faces")) {
        merge_faces = true;
    }

    if (!merge_faces) {
        return db::_3d::Surface_meshIO::load(sm, np);
    }

#if 1
    PolyhedronSPtr polyhedron = db::_3d::Surface_meshIO::load(sm, np);
    db::_3d::AbstractFile::mergeCoplanarFacets(polyhedron);
#else
    namespace PMP = CGAL::Polygon_mesh_processing;

    CGAL::Bbox_3 bbox = PMP::bbox(sm);
    const CGAL::FT diag_length = CGAL::approximate_sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                                        CGAL::square(bbox.ymax() - bbox.ymin()) +
                                                        CGAL::square(bbox.zmax() - bbox.zmin()));

    // Use shape detection to analyze the mesh
    std::vector<std::size_t> region_ids(num_faces(sm));
    boost::vector_property_map<Plane3> plane_map; // supporting planes of the regions detected

    const CGAL::FT cos_of_max_angle = 0.98;
    const CGAL::FT max_distance = 0.0001 * diag_length;

    // detect planar regions in the mesh
    // @todo growing should:
    // - use the .ini value of 'epsilon_coplanarity'
    // - stop if it merges faces with different weights
    // - give an error for adjacent coplanar faces that have different weights
    std::size_t nb_regions =
        PMP::region_growing_of_planes_on_faces(sm,
                                               CGAL::make_random_access_property_map(region_ids),
                                               CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                                .region_primitive_map(plane_map)
                                                                .maximum_distance(max_distance));

    static int region_dump_id = -1;
    utils::save_colored_mesh(sm, region_ids, "results/regions_" + std::to_string(++region_dump_id) + ".ply");

    // detect corner vertices on the boundary of planar regions
    std::vector<std::size_t> corner_ids(num_vertices(sm), -1); // corner status of vertices
    std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions

    std::size_t nb_corners =
        PMP::detect_corners_of_regions(sm,
                                      CGAL::make_random_access_property_map(region_ids),
                                      nb_regions,
                                      CGAL::make_random_access_property_map(corner_ids),
                                      CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle).
                                                        maximum_distance(max_distance).
                                                        edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

    CGAL_SS3_TRANSF_TRACE_CODE(for (face_descriptor f : faces(sm)))
    CGAL_SS3_TRANSF_TRACE("facet " << f << " is in region " << region_ids[f]);

    // the almost-coplanar merge is performed after the conversion to the Polyhedron
    // data structure because we want to be able to create faces that have holes,
    // which the CGAL::Surface_mesh class does not support
    std::map<edge_descriptor, EdgeWPtr> e2e;
    PolyhedronSPtr polyhedron = db::_3d::Surface_meshIO::load(sm, e2e, np);

    // merge the facets incident to an unconstrained edge (i.e., the edge is interior to a region)
    for (edge_descriptor e: edges(sm)) {
        if (ecm[e]) {
            continue;
        }

        EdgeSPtr edge = e2e[e].lock();
        if (!edge) {
            continue;
        }

        CGAL_SS3_TRANSF_TRACE("Merging facets " << edge->getFacetL()->getID() << " and " << edge->getFacetR()->getID());
        CGAL_assertion(sm.point(source(e, sm)) == *(edge->getVertexSrc()->getPoint()));
        CGAL_assertion(sm.point(target(e, sm)) == *(edge->getVertexDst()->getPoint()));

        // @fixme it seems like intermediate states are somewhat unsound during edge merging
        db::_3d::AbstractFile::mergeFacets(edge, polyhedron);
    }

    polyhedron->initializeAllIDs();

    db::_3d::AbstractFile::sanitize(polyhedron);
#endif

    CGAL_SS3_TRANSF_TRACE("Converted, " << polyhedron->facets().size() << " facets");

    return polyhedron;
}

} }

#endif /* ALGO_3D_POLYHEDRONTRANSFORMATION_H */
