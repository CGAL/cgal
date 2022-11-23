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

#ifndef CGAL_CARTESIAN_GRID_GEOMETRY_H
#define CGAL_CARTESIAN_GRID_GEOMETRY_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Grid_topology.h>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits>
class Cartesian_grid_geometry {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;

    typedef typename Grid_topology::Vertex_descriptor Vertex_descriptor;

public:
    Cartesian_grid_geometry(const Vector& offset, const Vector& spacing) : offset(offset), spacing(spacing) {}

    Point operator()(const Vertex_descriptor& v) const {
        return Point(v[0] * spacing[0], v[1] * spacing[1], v[2] * spacing[2]) + offset;
    }

private:
    Vector offset;
    Vector spacing;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_CARTESIAN_GRID_GEOMETRY_H
