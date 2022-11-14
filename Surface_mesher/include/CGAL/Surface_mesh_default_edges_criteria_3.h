// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_SURFACE_MESH_DEFAULT_EDGES_CRITERIA_3_H
#define CGAL_SURFACE_MESH_DEFAULT_EDGES_CRITERIA_3_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Surface_mesher/Surface_mesher_edges_level.h>

namespace CGAL {

template <class Tr, class Surface>
class Surface_mesh_default_edges_criteria_3
{
public:
  typedef typename Tr::Geom_traits GT;
  typedef typename GT::FT FT;
  typedef typename GT::Point_3 Point_3;

  typedef typename Tr::Edge Edge;
  typedef typename Tr::Vertex_handle Vertex_handle;

  Surface_mesh_default_edges_criteria_3(const FT radius_bound,
                                        const FT distance_bound,
                                        const Surface& surface,
                                        const bool assure_vertices_are_on_same_curve = true)
    : sq_distance_bound(distance_bound*distance_bound),
      sq_radius_bound(FT(4)*radius_bound*radius_bound),
      assure_vertices_are_on_same_curve(assure_vertices_are_on_same_curve),
      surface(surface)
  {
  }

  bool is_bad (const Edge& e,
               const Point_3& lineic_center) const
  {
    const Vertex_handle& va = e.first->vertex(e.second);
    const Vertex_handle& vb = e.first->vertex(e.third);

    if(assure_vertices_are_on_same_curve &&
       surface.vertices_not_on_same_curve(va, vb))
    {
      CGAL_MESHES_OUTPUT_STREAM << "e";
      return true;
    }

    typename GT::Compute_squared_distance_3 sq_distance =
      GT().compute_squared_distance_3_object();

    if(sq_radius_bound != FT(0) &&
       sq_distance(va->point(), vb->point()) > sq_radius_bound)
      return true;

    typename GT::Construct_midpoint_3 midpoint =
      GT().construct_midpoint_3_object();

    return sq_distance_bound != FT(0) &&
      sq_distance(lineic_center, midpoint(va->point(), vb->point())) >
      sq_distance_bound;
  }
private:
  FT sq_distance_bound;
  FT sq_radius_bound;
  bool assure_vertices_are_on_same_curve;
  const Surface& surface;
}; // end class Surface_mesh_default_edges_criteria_3

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_DEFAULT_EDGES_CRITERIA_3_H
