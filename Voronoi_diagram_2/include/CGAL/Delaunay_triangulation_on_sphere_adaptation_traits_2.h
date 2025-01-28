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

template <typename DToS2>
struct Delaunay_triangulation_on_sphere_adaptation_traits_2
{
public:
  typedef DToS2                                                      Delaunay_graph;
  typedef typename DToS2::Geom_traits                                Geom_traits;
  typedef typename Geom_traits::Point_on_sphere_2                    Point_2;
  typedef Point_2                                                    Site_2;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Point_accessor<Point_2,DToS2,Tag_true> Access_site_2;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::DToS2_Voronoi_point_2<DToS2>   Construct_Voronoi_point_2;

  typedef typename Delaunay_graph::Vertex_handle                     Delaunay_vertex_handle;
  typedef typename Delaunay_graph::Edge                              Delaunay_edge;
  typedef typename Delaunay_graph::Face_handle                       Delaunay_face_handle;

  typedef CGAL::Tag_false                                            Has_nearest_site_2;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Null_functor                   Nearest_site_2;

  Delaunay_triangulation_on_sphere_adaptation_traits_2(const Geom_traits& gt) : gt(gt) { }

  Access_site_2 access_site_2_object() const
  { return Access_site_2(gt); }

  Construct_Voronoi_point_2 construct_Voronoi_point_2_object() const
  { return Construct_Voronoi_point_2(gt); }

  Nearest_site_2 nearest_site_2_object() const
  { return Nearest_site_2(); }

private:
  const Geom_traits gt; // intentional copy
};

} //namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_ADAPTATION_TRAITS_2_H
