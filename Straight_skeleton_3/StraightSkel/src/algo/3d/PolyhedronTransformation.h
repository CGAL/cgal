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

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

class PolyhedronTransformation {
public:
    virtual ~PolyhedronTransformation();

    static void translate(PolyhedronSPtr polyhedron, Vector3SPtr t);
    static void scale(PolyhedronSPtr polyhedron, Vector3SPtr s);

    /**
     * returns the position of the vertex of a polyhedron, computed from the planes of
     * its incident faces.
     */
    static void resetPoint(VertexSPtr vertex);

    /**
     * returns the shifted position of the vertex of a polyhedron
     */
    static Point3SPtr shiftPoint(VertexSPtr vertex, CGAL::FT offset);

    /**
     * returns the shifted position of the facet of a polyhedron
     */
    static Plane3SPtr shiftPlane(FacetSPtr vertex, CGAL::FT offset);

    /**
     * Offsets the polyhedron `polyhedron`
     * Negative offset points to the interior of the polyhedron.
     */
    static void shiftFacetsInPlace(PolyhedronSPtr polyhedron,
                                   CGAL::FT offset,
                                   const bool recompute_positions = true);

    /**
     * Creates an offset polyhedron.
     * Negative offset points to the interior of the polyhedron.
     */
    static PolyhedronSPtr shiftFacets(PolyhedronSPtr polyhedron,
                                      CGAL::FT offset,
                                      const bool recompute_positions = true);


    /**
     * Normalize facet planes
    */
    static void normalizeFacetPlanes(PolyhedronSPtr polyhedron);

    /**
     * Normalize facet planes, ensuring parallel facets receive the same plane coefficients.
    */
    static void harmonizeFacetPlanes(PolyhedronSPtr polyhedron);

    /**
     * To check for parallel planes is not enough.
     */
    static bool hasParallelPlanes(PolyhedronSPtr polyhedron);
    static bool doAll3PlanesIntersect(PolyhedronSPtr polyhedron);
    static void randMovePoints(PolyhedronSPtr polyhedron);
    static PolyhedronSPtr perturb(PolyhedronSPtr polyhedron);

    static Point3SPtr boundingBoxMin(PolyhedronSPtr polyhedron);
    static Point3SPtr boundingBoxMax(PolyhedronSPtr polyhedron);

    static void translateNscale(PolyhedronSPtr polyhedron,
            Point3SPtr p_box_min, Point3SPtr p_box_max);

    static bool isInsideBox(PolyhedronSPtr polyhedron,
            Point3SPtr p_box_min, Point3SPtr p_box_max);

protected:
    static Vector3SPtr randVec(double min, double max);
    PolyhedronTransformation();
};

} }

#endif /* ALGO_3D_POLYHEDRONTRANSFORMATION_H */
