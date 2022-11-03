// Copyright (c) 2022 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_OCTREE_GEOMETRY_H
#define CGAL_OCTREE_GEOMETRY_H

#include <CGAL/Isosurfacing_3/internal/Octree_topology.h>
#include <CGAL/Octree_wrapper.h>
#include <CGAL/license/Isosurfacing_3.h>

#include <memory>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits>
class Octree_geometry {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::Point_3 Point;

    typedef std::shared_ptr<Octree_wrapper<Geom_traits>> Octree;

    typedef typename Octree_topology<Geom_traits>::Vertex_descriptor Vertex_descriptor;

public:
    Octree_geometry(const Octree& octree) : octree(octree) {}

    Point operator()(const Vertex_descriptor& v) const {
        return octree->point(v);
    }

private:
    const Octree octree;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_OCTREE_GEOMETRY_H
