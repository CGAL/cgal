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
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_APOLLONIUS_GRAPH_NEAREST_SITE_2_H
#define CGAL_VORONOI_DIAGRAM_2_APOLLONIUS_GRAPH_NEAREST_SITE_2_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Apollonius_graph_2/Orientation_2.h>

#include <boost/variant.hpp>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================

template<class DG>
class Apollonius_graph_nearest_site_2
{
 public:
  typedef DG                                          Delaunay_graph;
  typedef typename Delaunay_graph::Point_2            Point_2;

 private:
  typedef Triangulation_cw_ccw_2                      CW_CCW_2;

  typedef typename Delaunay_graph::Geom_traits        Geom_traits;
  typedef typename Delaunay_graph::Site_2             Site_2;
  typedef typename Delaunay_graph::Vertex_handle      Vertex_handle;
  typedef typename Delaunay_graph::Face_handle        Face_handle;
  typedef typename Delaunay_graph::Edge               Edge;
  typedef typename Delaunay_graph::Face_circulator    Face_circulator;
  typedef typename Delaunay_graph::Edge_circulator    Edge_circulator;

public:
  typedef boost::variant<Vertex_handle,Edge,Face_handle> result_type;

  result_type operator()(const Delaunay_graph& dg, const Point_2& p) const {
    CGAL_precondition( dg.dimension() >= 0 );

    typename Geom_traits::Oriented_side_of_bisector_2 side_of_bisector =
      dg.geom_traits().oriented_side_of_bisector_2_object();

    Vertex_handle v = dg.nearest_neighbor(p);
    if ( dg.dimension() == 0 ) {
      return v;
    }

    if ( dg.dimension() == 1 ) {
      Edge e = *dg.finite_edges_begin();
      Vertex_handle v1 = e.first->vertex(CW_CCW_2::ccw(e.second));
      Vertex_handle v2 = e.first->vertex(CW_CCW_2::cw(e.second) );

      Oriented_side os = side_of_bisector(v1->site(), v2->site(), p);
      
      if ( os == ON_ORIENTED_BOUNDARY ) {
	return e;
      } else {
	return v;
      }
    }

    CGAL_assertion( dg.dimension() == 2 );

    Face_circulator fc_start = dg.incident_faces(v);
    Face_circulator fc = fc_start;

    // first check if the point lies on a Voronoi vertex
    do {
      int index = fc->index(v);
      Vertex_handle v1 = fc->vertex(CW_CCW_2::ccw(index));
      Vertex_handle v2 = fc->vertex(CW_CCW_2::cw(index) );

      Oriented_side os1 = ON_POSITIVE_SIDE, os2 = ON_POSITIVE_SIDE;

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	os1 = side_of_bisector(v->site(), v1->site(), p);
      }
      if ( !dg.is_infinite(v2) ) {
	os2 = side_of_bisector(v->site(), v2->site(), p);
      }

      CGAL_assertion( os1 != ON_NEGATIVE_SIDE );
      CGAL_assertion( os2 != ON_NEGATIVE_SIDE );

      if ( os1 == ON_ORIENTED_BOUNDARY && os2 == ON_ORIENTED_BOUNDARY ) {
	return Face_handle(fc);
      }

      ++fc;
    } while ( fc != fc_start );

    // now check if the point lies on a Voronoi edge
    Edge_circulator ec_start = dg.incident_edges(v);
    Edge_circulator ec = ec_start;
    do {
      Face_handle f = ec->first;
      int idx = ec->second;
      CGAL_assertion( f->vertex(CW_CCW_2::cw(idx)) == v );
      Vertex_handle v1 = f->vertex(CW_CCW_2::ccw(idx));

      Oriented_side os1 = ON_POSITIVE_SIDE;

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	os1 = side_of_bisector(v->site(), v1->site(), p);
      }

      CGAL_assertion( os1 != ON_NEGATIVE_SIDE );

      if ( os1 != ON_ORIENTED_BOUNDARY ) {
	++ec;
	continue;
      }

      // now find the correct part of the bisector on which the query
      // point lies; this part is essential when the two sites have
      // more than one Voronoi edge between them
      Vertex_handle v2 = f->vertex(idx);
      Vertex_handle v3 = dg.tds().mirror_vertex(f, idx);

      bool is_inf2 = dg.is_infinite(v2);
      bool is_inf3 = dg.is_infinite(v3);

      typename Geom_traits::Orientation_2 vv_orientation =
	dg.geom_traits().orientation_2_object();

      Orientation o2 = LEFT_TURN, o3 = RIGHT_TURN;
      if ( is_inf2 && is_inf3 ) { return *ec; }

      Site_2 sp(p, 0);

      if ( !is_inf2 && is_inf3 ) {
	Orientation vo2 =
	  vv_orientation(v->site(), v2->site(), v1->site(),
			 v->site(), v1->site());
	Orientation op = vv_orientation(v->site(), v1->site(), sp);

	if ( vo2 != LEFT_TURN && op == LEFT_TURN ) { return *ec; }

	o2 = vv_orientation(v->site(), v2->site(), v1->site(),
			    v->site().point(), sp);
	CGAL_assertion( o2 != COLLINEAR );
	if ( o2 == RIGHT_TURN ) { return *ec; }
      }

      if ( is_inf2 && !is_inf3 ) {
	Orientation vo3 =
	  vv_orientation(v->site(), v1->site(), v3->site(),
			 v->site(), v1->site());
	Orientation op = vv_orientation(v->site(), v1->site(), sp);

	if ( vo3 != RIGHT_TURN && op == RIGHT_TURN ) { return *ec; }

	o3 = vv_orientation(v->site(), v1->site(), v3->site(),
			    v->site().point(), sp);
	CGAL_assertion( o3 != COLLINEAR );
	if ( o3 == LEFT_TURN ) { return *ec; }
      }

      if ( !is_inf2 && !is_inf3 ) {
	Orientation vo2 =
	  vv_orientation(v->site(), v2->site(), v1->site(),
			 v->site(), v1->site());
	Orientation vo3 =
	  vv_orientation(v->site(), v1->site(), v3->site(),
			 v->site(), v1->site());

	if ( (vo3 == LEFT_TURN && vo2 != RIGHT_TURN) ||
	     (vo3 != LEFT_TURN && vo2 == RIGHT_TURN) ) {
	  o2 = vv_orientation(v->site(), v2->site(), v1->site(),
			      v->site(), sp);
	  o3 = vv_orientation(v->site(), v1->site(), v3->site(),
			      v->site(), sp);
	  CGAL_assertion( o2 != COLLINEAR && o3 != COLLINEAR );
	  if ( o2 == RIGHT_TURN && o3 == LEFT_TURN ) { return *ec; }
	} else {
	  CGAL_assertion( vo2 == RIGHT_TURN && vo3 == LEFT_TURN );

	  Orientation op = vv_orientation(v->site(), v1->site().point(), sp);

	  if ( op == COLLINEAR ) { return *ec; }
	  else if ( op == LEFT_TURN ) {
	    o3 = vv_orientation(v->site(), v1->site(), v3->site(),
				v->site(), sp);
	    CGAL_assertion( o3 != COLLINEAR );
	    if ( o3 == LEFT_TURN ) { return *ec; }
	  } else {
	    o2 = vv_orientation(v->site(), v2->site(), v1->site(),
				v->site(), sp);
	    CGAL_assertion( o2 != COLLINEAR );
	    if ( o2 == RIGHT_TURN ) { return *ec; }
	  }
	}
      }
      ++ec;
    } while ( ec != ec_start );

    // the point lies in a Voronoi face
    return v;
  }
};


//=========================================================================
//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_APOLLONIUS_GRAPH_NEAREST_SITE_2_H
