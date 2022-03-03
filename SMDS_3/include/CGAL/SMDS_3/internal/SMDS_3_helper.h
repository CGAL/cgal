// Copyright (c) 2021 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_INTERNAL_SMDS_3_HELPERS_H
#define CGAL_INTERNAL_SMDS_3_HELPERS_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/enum.h>

#include <vector>
#include <array>

namespace CGAL {
namespace SMDS_3 {
namespace internal {

  template<typename Triangulation>
  bool is_convex(const Triangulation& tr)
  {
    typedef typename Triangulation::Vertex_handle        Vertex_handle;
    typedef typename Triangulation::Cell_handle          Cell_handle;
    typedef typename Triangulation::Geom_traits::Point_3 Point_3;
    typename Triangulation::Geom_traits::Construct_point_3 cp
      = tr.geom_traits().construct_point_3_object();
    typename Triangulation::Geom_traits::Orientation_3 orientation =
      tr.geom_traits().orientation_3_object();

    std::vector<Cell_handle> infcells;
    tr.incident_cells(tr.infinite_vertex(), std::back_inserter(infcells));
    for (Cell_handle c : infcells)
    {
      const Cell_handle neigh = c->neighbor(c->index(tr.infinite_vertex()));
      const int i = neigh->index(c);

      const std::array<Point_3, 3> pfacet = { cp(neigh->vertex((i + 1) % 4)->point()),
                                              cp(neigh->vertex((i + 2) % 4)->point()),
                                              cp(neigh->vertex((i + 3) % 4)->point())};
      const CGAL::Orientation o = orientation(
          pfacet[0], pfacet[1], pfacet[2], cp(neigh->vertex(i)->point()));

      for (Vertex_handle v : tr.finite_vertex_handles())
      {
        if (c->has_vertex(v))
          continue;
        if (o != orientation(pfacet[0], pfacet[1], pfacet[2],
                             cp(neigh->vertex(i)->point())))
          return false;
      }
    }

    return true;
  }

} // end namespace internal
} // end namespace SMDS_3
} // end namespace CGAL

#endif // CGAL_INTERNAL_SMDS_3_HELPERS_H
