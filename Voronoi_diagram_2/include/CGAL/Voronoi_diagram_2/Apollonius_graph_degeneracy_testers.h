// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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

#ifndef CGAL_VORONOI_DIAGRAM_2_APOLLONIUS_GRAPH_DEGENERACY_TESTERS_H
#define CGAL_VORONOI_DIAGRAM_2_APOLLONIUS_GRAPH_DEGENERACY_TESTERS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_base_2.h>
#include <CGAL/Voronoi_diagram_2/Identity_rejectors.h>
#include <CGAL/Triangulation_utils_2.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================

template<class DG>
class Apollonius_graph_edge_tester_2
  : public Rejector_base
{
  // tests whether a dual edge has zero length
 public:
  typedef DG                                       Delaunay_graph;

  typedef typename DG::Edge                        Edge;
  typedef typename DG::Face_handle                 Face_handle;
  typedef typename DG::Edge_circulator             Edge_circulator;
  typedef typename DG::All_edges_iterator          All_edges_iterator;
  typedef typename DG::Finite_edges_iterator       Finite_edges_iterator;

  typedef bool           result_type;

 private:
  typedef Triangulation_cw_ccw_2                     CW_CCW_2;

  typedef typename Delaunay_graph::Geom_traits       Geom_traits;

  typedef typename Delaunay_graph::Vertex_handle     Vertex_handle;

  typedef typename Geom_traits::Site_2               Site_2;

 public:
  bool operator()(const Delaunay_graph& dual,
                  const Face_handle& f, int i) const
  {
    if ( dual.dimension() == 1 ) { return false; }

    if ( dual.is_infinite(f, i) ) { return false; }

    Vertex_handle v3 = f->vertex(     i  );
    Vertex_handle v4 = dual.data_structure().mirror_vertex(f, i);

    if ( dual.is_infinite(v3) || dual.is_infinite(v4) ) {
      return false;
    }

    Vertex_handle v1 = f->vertex( CW_CCW_2::ccw(i) );
    Vertex_handle v2 = f->vertex( CW_CCW_2::cw(i) );

    Site_2 s1 = v1->site();
    Site_2 s2 = v2->site();
    Site_2 s3 = v3->site();
    Site_2 s4 = v4->site();
    return dual.geom_traits().is_degenerate_edge_2_object()(s1,s2,s3,s4);
  }

  bool operator()(const Delaunay_graph& dual, const Edge& e) const {
    return operator()(dual, e.first, e.second);
  }

  bool operator()(const Delaunay_graph& dual,
                  const All_edges_iterator& eit) const {
    return operator()(dual, *eit);
  }

  bool operator()(const Delaunay_graph& dual,
                  const Finite_edges_iterator& eit) const {
    return operator()(dual, *eit);
  }

  bool operator()(const Delaunay_graph& dual,
                  const Edge_circulator& ec) const {
    return operator()(dual, *ec);
  }
};


//=========================================================================
//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_APOLLONIUS_GRAPH_DEGENERACY_TESTERS_H
