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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_DELAUNAY_TRIANGULATION_ADAPTATION_TRAITS_2_H
#define CGAL_DELAUNAY_TRIANGULATION_ADAPTATION_TRAITS_2_H 1

#include <CGAL/license/Voronoi_diagram_2.h>

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_functors.h>
#include <CGAL/Voronoi_diagram_2/Construct_dual_points.h>
#include <CGAL/Voronoi_diagram_2/Site_accessors.h>

namespace CGAL {

template <typename DTOS>
struct Delaunay_triangulation_on_sphere_adaptation_traits_2
{
public:
  typedef DTOS                                                       Delaunay_graph;
  typedef typename Delaunay_graph::Geom_traits                       Geom_traits;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::DToS2_Point_accessor<DTOS>     Access_site_2;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::DToS2_Voronoi_point_2<DTOS>    Construct_Voronoi_point_2;

  typedef typename Delaunay_graph::Vertex_handle                     Delaunay_vertex_handle;
  typedef typename Delaunay_graph::Edge                              Delaunay_edge;
  typedef typename Delaunay_graph::Face_handle                       Delaunay_face_handle;

  typedef CGAL::Tag_false                                            Has_nearest_site_2;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Null_functor                   Nearest_site_2;

  Delaunay_triangulation_on_sphere_adaptation_traits_2(const DTOS& dtos) : dtos(dtos) { }

  Access_site_2 access_site_2_object() const
  { return Access_site_2(dtos); }

  Construct_Voronoi_point_2 construct_Voronoi_point_2_object() const
  { return Construct_Voronoi_point_2(dtos); }

  Nearest_site_2 nearest_site_2_object() const
  { return Nearest_site_2(); }

  typedef typename Geom_traits::Point_on_sphere_2                    Point_2;
  typedef Point_2                                                    Site_2;

private:
  const DTOS& dtos;
};

} //namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_ADAPTATION_TRAITS_2_H
