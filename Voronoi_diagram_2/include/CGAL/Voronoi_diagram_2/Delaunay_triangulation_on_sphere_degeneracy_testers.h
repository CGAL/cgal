// Copyright (c) 2021 GeometryFactory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>
//                 Mael Rouxel-Labbé

#ifndef CGAL_VORONOI_DIAGRAM_2_DELAUNAY_TRIANGULATION_ON_SPHERE_DEGENERACY_TESTERS_H
#define CGAL_VORONOI_DIAGRAM_2_DELAUNAY_TRIANGULATION_ON_SPHERE_DEGENERACY_TESTERS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_base_2.h>
#include <CGAL/Voronoi_diagram_2/Identity_rejectors.h>

namespace CGAL {
namespace VoronoiDiagram_2 {
namespace Internal {

//=========================================================================
//=========================================================================

// tests whether a dual edge has zero length
template <class DG>
class Delaunay_triangulation_on_sphere_edge_tester_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Rejector_base
{
public:
  typedef DG                                             Delaunay_graph;

  typedef typename Delaunay_graph::Edge                  Edge;
  typedef typename Delaunay_graph::Face_handle           Face_handle;
  typedef typename Delaunay_graph::Edge_circulator       Edge_circulator;
  typedef typename Delaunay_graph::All_edges_iterator    All_edges_iterator;
  typedef typename Delaunay_graph::Finite_edges_iterator Finite_edges_iterator;
  typedef bool                                           result_type;

 private:
  typedef typename Delaunay_graph::Geom_traits           Geom_traits;
  typedef typename Delaunay_graph::Vertex_handle         Vertex_handle;
  typedef typename Delaunay_graph::Point                 Point;

 public:
  bool operator()(const Delaunay_graph& dual,
                  const Face_handle f, int i) const
  {
    if(dual.dimension() == 1)
      return false;

    const Vertex_handle v3 = f->vertex(i);
    const Vertex_handle v4 = dual.tds().mirror_vertex(f, i);

    const Vertex_handle v1 = f->vertex(dual.ccw(i));
    const Vertex_handle v2 = f->vertex(dual.cw(i));

    const Point& p1 = v1->point();
    const Point& p2 = v2->point();
    const Point& p3 = v3->point();
    const Point& p4 = v4->point();
    const Oriented_side os = dual.geom_traits().side_of_oriented_circle_on_sphere_2_object()(p1,p2,p3,p4);

    return (os == ON_ORIENTED_BOUNDARY);
  }

  bool operator()(const Delaunay_graph& dual, const Edge& e) const
  {
    return operator()(dual, e.first, e.second);
  }

  template <typename EdgeIterator>
  bool operator()(const Delaunay_graph& dual,
                  const EdgeIterator eit) const
  {
    return operator()(dual, *eit);
  }
};

//=========================================================================
//=========================================================================

} // namespace Internal
} // namespace VoronoiDiagram_2
} // namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_DELAUNAY_TRIANGULATION_ON_SPHERE_DEGENERACY_TESTERS_H
