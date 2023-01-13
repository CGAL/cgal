// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_INTERNAL_OCTREE_GEOMETRY_H
#define CGAL_ISOSURFACING_3_INTERNAL_OCTREE_GEOMETRY_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Octree_topology.h>
#include <CGAL/Isosurfacing_3/internal/Octree_wrapper.h>

#include <memory>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <typename Octree>
class Octree_geometry
{
  using Vertex_descriptor = typename Octree_topology<Octree>::Vertex_descriptor;

public:
  Octree_geometry(const Octree& octree)
    : octree(octree)
  { }

  decltype(auto) /*Point_3*/ operator()(const Vertex_descriptor& v) const
  {
    return octree.point(v);
  }

private:
  const Octree& octree;
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_OCTREE_GEOMETRY_H
