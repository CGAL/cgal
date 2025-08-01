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
 * @file   algo/3d/SelfIntersection.h
 * @author Gernot Walzl
 * @date   2012-07-18
 */

#ifndef ALGO_3D_SELFINTERSECTION_H
#define ALGO_3D_SELFINTERSECTION_H

#include "data/3d/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;

class SelfIntersection {
public:
    virtual ~SelfIntersection();

    static bool doEdgesShareAVertex(EdgeSPtr edge1, EdgeSPtr edge2, bool handle_deg1_as_ray);
    static bool doEdgesIntersect(FacetSPtr facet, EdgeSPtr edge1, EdgeSPtr edge2, bool handle_deg1_as_ray);
    static bool isSelfIntersectingFacet(FacetSPtr facet);
    static unsigned int hasSelfIntersectingFacets(PolyhedronSPtr polyhedron);

    static bool isInsideWithRayShooting(const Point3& point,
                                        FacetSPtr facet,
                                        const bool handle_deg1_as_ray);
    static bool isInsideWithRayShootingV2(Point3SPtr point,
                                          FacetSPtr facet);
    static bool isEdgeInsideFacet(FacetSPtr facet, EdgeSPtr edge, bool handle_deg1_as_ray);
    static bool hasSelfIntersectingSurface(PolyhedronSPtr polyhedron);

protected:
    SelfIntersection();
};

} }

#endif /* ALGO_3D_SELFINTERSECTION_H */
