// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_SEGMENT_DELAUNAY_GRAPH_DEGENERACY_TESTERS_H
#define CGAL_VORONOI_DIAGRAM_2_SEGMENT_DELAUNAY_GRAPH_DEGENERACY_TESTERS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_base_2.h>
#include <CGAL/Voronoi_diagram_2/Identity_rejectors.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================

template<class DG>
class Segment_Delaunay_graph_edge_tester_2
  : public Rejector_base
{
  // tests whether a dual edge has zero length
 public:
  typedef DG                                           Delaunay_graph;

  typedef typename Delaunay_graph::Edge                Edge;
  typedef typename Delaunay_graph::Face_handle         Face_handle;
  typedef typename Delaunay_graph::Edge_circulator     Edge_circulator;
  typedef typename Delaunay_graph::All_edges_iterator  All_edges_iterator;

  typedef typename Delaunay_graph::Finite_edges_iterator
  Finite_edges_iterator;

  typedef bool           result_type;

 private:
  typedef typename Delaunay_graph::Geom_traits         Geom_traits;
  typedef typename Delaunay_graph::Vertex_handle       Vertex_handle;
  typedef typename Delaunay_graph::Site_2              Site_2;

  typedef typename Geom_traits::Equal_2                Equal_2;
  typedef typename Geom_traits::Orientation_2          Orientation_2;

 private:
  bool is_degenerate_infinite_edge(const Delaunay_graph& dual,
				   const Face_handle& f, int i) const
  {
    CGAL_precondition( dual.is_infinite(f, i) );

    Vertex_handle v = f->vertex( dual.ccw(i) );
    Vertex_handle v_inf = f->vertex( dual.cw(i) );

    if ( dual.is_infinite(v) ) {
      std::swap(v, v_inf);
    }

    if ( v->storage_site().is_segment() ) { return false; }

    Vertex_handle vv[2];

    vv[0] = f->vertex(i);
    vv[1] = dual.data_structure().mirror_vertex(f, i);

    if ( vv[0] == vv[1] ) { return false; }

    Site_2 s_end[2];
    for (int i = 0; i < 2; i++) {
      if ( vv[i]->storage_site().is_point() ) { return false; }

      Equal_2 are_equal = dual.geom_traits().equal_2_object();
      if ( !are_equal(v->site(),vv[i]->site().source_site()) &&
	   !are_equal(v->site(),vv[i]->site().target_site()) ) {
	return false;
      }
      if ( are_equal(v->site(),vv[i]->site().source_site()) ) {
	s_end[i] = vv[i]->site().target_site();
      } else {
	s_end[i] = vv[i]->site().source_site();
      }
    }

    Orientation_2 orientation = dual.geom_traits().orientation_2_object();

    return orientation(s_end[0], s_end[1], v->site()) == COLLINEAR;
  }

 public:
  bool operator()(const Delaunay_graph& dual,
		  const Face_handle& f, int i) const
  {
    if ( dual.dimension() == 1 ) { return false; }

    if ( dual.is_infinite(f, i) ) {
      return is_degenerate_infinite_edge(dual, f, i);
    }

    Vertex_handle v3 = f->vertex(i);
    Vertex_handle v4 = dual.data_structure().mirror_vertex(f, i);

    if ( dual.is_infinite(v3) || dual.is_infinite(v4) ) {
      return false;
    }

    Vertex_handle v1 = f->vertex( dual.ccw(i) );
    Vertex_handle v2 = f->vertex( dual.cw(i) );

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

template<class DG>
class Segment_Delaunay_graph_face_tester_2
  : public Rejector_base
{
  // tests whether a face has zero area
 public:
  typedef DG                                       Delaunay_graph;
  typedef typename Delaunay_graph::Vertex_handle   Vertex_handle;

  typedef bool           result_type;

 private:
  typedef Segment_Delaunay_graph_edge_tester_2<Delaunay_graph>  Edge_tester;

  typedef typename Delaunay_graph::Geom_traits      Geom_traits;
  typedef typename Delaunay_graph::Edge             Edge;
  typedef typename Delaunay_graph::Edge_circulator  Edge_circulator;
  typedef typename Delaunay_graph::Face_handle      Face_handle;
  typedef typename Delaunay_graph::size_type        size_type;
  typedef typename Delaunay_graph::Site_2           Site_2;

  typedef typename Geom_traits::Equal_2             Equal_2;
  typedef typename Geom_traits::Orientation_2       Orientation_2;

 public:
  bool operator()(const Delaunay_graph& dual, const Vertex_handle& v) const
  {
    if ( dual.dimension() < 2 ) { return false; }

    if ( dual.is_infinite(v) ) { return false; }

    // THIS TEST NEEDS TO USE GEOMETRY; I CANNOT DO IT IN AN ENTIRELY
    // COMBINATORIAL MANNER

    // SEGMENT SPECIFIC TEST
    if ( v->site().is_segment() ) { return false; }

    // THIS WORKS ONLY FOR SEGMENTS (OR MAYBE NOT...)
    Edge_circulator ec_start(v);
    Edge_circulator ec = ec_start;
    size_type deg = 0;       // vertex degree
    size_type n_degen = 0;   // number of degenerate/non-infinite edges
    size_type n_inf = 0;     // number of infinite edges
    // number of non-degenerate/non-infinite edges
    size_type n_non_degen = 0;
      
    Edge e[2];
    Edge_tester e_tester;
    do {
      if ( e_tester(dual, ec) ) { ++n_degen; }
      else if ( dual.is_infinite(ec) ) { ++n_inf; }
      else { 
	if ( !dual.is_infinite(ec) ) {
	  if ( n_non_degen < 2 ) {
	    e[n_non_degen] = *ec;
	  }
	  n_non_degen++;
	}
      }
      deg++;
      ++ec;
    } while ( ec != ec_start );

    if ( deg == n_degen ) { return true; }
    if ( n_non_degen != 2 ) { return false; }

    Vertex_handle vv[2];
    Site_2 s_end[2];
    for (int i = 0; i < 2; i++) {
      CGAL_assertion( !dual.is_infinite(e[i]) );
      CGAL_expensive_assertion( !e_tester(dual, e[i]) );

      Vertex_handle v1 = e[i].first->vertex( dual.ccw(e[i].second) );
      Vertex_handle v2 = e[i].first->vertex( dual.cw(e[i].second) );
      vv[i] = (v1 == v) ? v2 : v1;

      CGAL_assertion( v == v1 || v == v2 );

      if ( vv[i]->storage_site().is_point() ) { return false; }

      Equal_2 are_equal = dual.geom_traits().equal_2_object();
      if ( !are_equal(v->site(),vv[i]->site().source_site()) &&
	   !are_equal(v->site(),vv[i]->site().target_site()) ) {
	return false;
      }
      if ( are_equal(v->site(),vv[i]->site().source_site()) ) {
	s_end[i] = vv[i]->site().target_site();
      } else {
	s_end[i] = vv[i]->site().source_site();
      }
    }

    Orientation_2 orientation = dual.geom_traits().orientation_2_object();

    return orientation(s_end[0], s_end[1], v->site()) == COLLINEAR;
  }

};

//=========================================================================
//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_SEGMENT_DELAUNAY_GRAPH_DEGENERACY_TESTERS_H
