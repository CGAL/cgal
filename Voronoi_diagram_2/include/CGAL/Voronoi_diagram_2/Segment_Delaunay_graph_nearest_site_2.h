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

#ifndef CGAL_VORONOI_DIAGRAM_2_SEGMENT_DELAUNAY_GRAPH_NEAREST_SITE_2_H
#define CGAL_VORONOI_DIAGRAM_2_SEGMENT_DELAUNAY_GRAPH_NEAREST_SITE_2_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Triangulation_utils_2.h>

#include <boost/variant.hpp>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================

template<class DG>
class Segment_Delaunay_graph_nearest_site_2
{
public:
  typedef DG                                          Delaunay_graph;
  typedef typename Delaunay_graph::Point_2            Point_2;

private:
  typedef Triangulation_cw_ccw_2                      CW_CCW_2;

  typedef typename Delaunay_graph::Site_2             Site_2;
  typedef typename Delaunay_graph::Vertex_handle      Vertex_handle;
  typedef typename Delaunay_graph::Face_handle        Face_handle;
  typedef typename Delaunay_graph::Face_circulator    Face_circulator;
  typedef typename Delaunay_graph::Edge_circulator    Edge_circulator;
  typedef typename Delaunay_graph::Edge               Edge;
  typedef typename Delaunay_graph::Geom_traits        Geom_traits;

  typedef typename Geom_traits::Equal_2               Equal_2;
  typedef typename Geom_traits::Orientation_2         Orientation_2;
  typedef typename Geom_traits::Oriented_side_2       Oriented_side_2;
  typedef typename Geom_traits::Oriented_side_of_bisector_2
  Oriented_side_of_bisector_2;

private:
  static Site_2 opposite(const Site_2& t) {
    CGAL_precondition( t.is_segment() );
    if ( t.is_input() ) {
      return Site_2::construct_site_2(t.target_of_supporting_site(),
				      t.source_of_supporting_site());
    } else {
      if ( t.is_input(0) ) {
	return Site_2::construct_site_2(t.target_of_supporting_site(),
					t.source_of_supporting_site(),
					t.target_of_crossing_site(1),
					t.source_of_crossing_site(1),
					false);
      } else if ( t.is_input(1) ) {
	return Site_2::construct_site_2(t.target_of_supporting_site(),
					t.source_of_supporting_site(),
					t.target_of_crossing_site(0),
					t.source_of_crossing_site(0),
					true);
      } else {
	return Site_2::construct_site_2(t.target_of_supporting_site(),
					t.source_of_supporting_site(),
					t.target_of_crossing_site(1),
					t.source_of_crossing_site(1),
					t.target_of_crossing_site(0),
					t.source_of_crossing_site(0));
      }
    }
  }

  Site_2 orient_segment(const Site_2& s, const Site_2& p,
			const Orientation_2& orientation) const
  {
    CGAL_precondition( s.is_segment() );
    CGAL_precondition( p.is_point() );

    Orientation o = orientation(s.source_site(), s.target_site(), p);
    CGAL_assertion( o != COLLINEAR );
    if ( o == LEFT_TURN ) { return s; }
    return opposite(s);
  }

public:
  typedef boost::variant<Vertex_handle,Edge,Face_handle>  result_type;

  result_type operator()(const Delaunay_graph& dg, const Point_2& p) const {
    CGAL_precondition( dg.dimension() >= 0 );

    Oriented_side_of_bisector_2 side_of_bisector =
      dg.geom_traits().oriented_side_of_bisector_2_object();

    Equal_2 equal = dg.geom_traits().equal_2_object();

    Site_2 sp = Site_2::construct_site_2(p);

    Vertex_handle v = dg.nearest_neighbor(p);
    if ( dg.dimension() == 0 ) {
      return v;
    }

    if ( dg.dimension() == 1 ) {
      Edge e = *dg.finite_edges_begin();
      Vertex_handle v1 = e.first->vertex(CW_CCW_2::ccw(e.second));
      Vertex_handle v2 = e.first->vertex(CW_CCW_2::cw(e.second) );

      Oriented_side os = side_of_bisector(v1->site(), v2->site(), sp);
      
      if ( os == ON_ORIENTED_BOUNDARY ) {
	return e;
      } else {
	return v;
      }
    }

    CGAL_assertion( dg.dimension() == 2 );

    Site_2 sv = v->site();

    Face_circulator fc_start = dg.incident_faces(v);
    Face_circulator fc = fc_start;

    // first check if the point lies on a Voronoi vertex
    do {
      int index = fc->index(v);
      Vertex_handle v1 = fc->vertex(CW_CCW_2::ccw(index));
      Vertex_handle v2 = fc->vertex(CW_CCW_2::cw(index) );

      Oriented_side os1 = ON_POSITIVE_SIDE, os2 = ON_POSITIVE_SIDE;

      // check if the query point is identical to an endpoint of a
      // segment that has a Voronoi face with zero area.
      if ( !dg.is_infinite(v1) && !dg.is_infinite(v2) ) {
	Site_2 s1 = v1->site();
	Site_2 s2 = v2->site();
	if ( sv.is_point() && equal(sv, sp) &&
	     s1.is_segment() && s2.is_segment() ) {
	  bool b1 = equal(sv, s1.source_site()) || equal(sv, s1.target_site());
	  bool b2 = equal(sv, s2.source_site()) || equal(sv, s2.target_site());

	  if ( b1 && b2 ) { return Face_handle(fc); }
	}
      }

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	os1 = side_of_bisector(sv, v1->site(), sp);
      }
      if ( !dg.is_infinite(v2) ) {
	os2 = side_of_bisector(sv, v2->site(), sp);
      }

      CGAL_assertion( os1 != ON_NEGATIVE_SIDE );
      CGAL_assertion( os2 != ON_NEGATIVE_SIDE );

      if ( os1 == ON_ORIENTED_BOUNDARY && os2 == ON_ORIENTED_BOUNDARY ) {
	return fc;
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

      // check if the query point is lies on the bisector between a
      // segment and its endpoint
      if ( !dg.is_infinite(v1) ) {
	Site_2 s1 = v1->site();
	if ( sv.is_point() &&  s1.is_segment() ) {
	  bool b = equal(sv, s1.source_site()) || equal(sv, s1.target_site());

	  if ( b && equal(sv, sp) ) { return *ec; }

	  Oriented_side os = side_of_bisector(sv, s1, sp);

	  if ( b && os == ON_ORIENTED_BOUNDARY ) { return *ec; }
	}
	if ( sv.is_segment() && s1.is_point() ) {
	  bool b = equal(sv.source_site(), s1) || equal(sv.target_site(), s1);

	  if ( b && equal(s1, sp) ) { return *ec; }

	  Oriented_side os = side_of_bisector(sv, s1, sp);

	  if ( b && os == ON_ORIENTED_BOUNDARY ) { return *ec; }
	}
      }

      Oriented_side os1 = ON_POSITIVE_SIDE;

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	os1 = side_of_bisector(sv, v1->site(), sp);
      }

      CGAL_assertion( os1 != ON_NEGATIVE_SIDE );

      if ( os1 != ON_ORIENTED_BOUNDARY ) {
	++ec;
	continue;
      }

      // now find the correct part of the bisector on which the query
      // point lies; this part is essential when the two sites have
      // more than one Voronoi edge between them
      CGAL_assertion( !dg.is_infinite(v1) );
      Site_2 s1 = v1->site();

      // if both v and v1 are points there is no need to check for
      // other Voronoi edges; there is only one.
      if ( sv.is_point() && s1.is_point() ) {
	return *ec;
      }

      // the vertices v2 and v3 can only be infinite in two cases:
      // 1. when both v and v1 are points
      // 2. when v and v1 are a point and a segment and the point is
      //    an endpoint of the segment
      // none of the two cases can happen at this point in the code,
      // since we already have treated these two cases above.

      Vertex_handle v2 = f->vertex(idx);
      Vertex_handle v3 = dg.tds().mirror_vertex(f, idx);

      CGAL_assertion( !dg.is_infinite(v2) && !dg.is_infinite(v3) );


      Oriented_side_2 vv_oriented_side =
	dg.geom_traits().oriented_side_2_object();

      Orientation_2 orientation = dg.geom_traits().orientation_2_object();

      Oriented_side o2 = ON_POSITIVE_SIDE, o3 = ON_NEGATIVE_SIDE;

      Site_2 sp = Site_2::construct_site_2(p);

      if ( sv.is_segment() ) {
	Site_2 svo = orient_segment(sv, sp, orientation);
	o2 = vv_oriented_side(sv, v2->site(), s1, svo, sp);
	o3 = vv_oriented_side(sv, s1, v3->site(), svo, sp);
	CGAL_assertion( o2 != ON_ORIENTED_BOUNDARY );
	CGAL_assertion( o3 != ON_ORIENTED_BOUNDARY );
	if ( o2 == ON_NEGATIVE_SIDE && o3 == ON_POSITIVE_SIDE ) {
	  return *ec;
	}
      } else {
	CGAL_assertion( s1.is_segment() );
	Site_2 s1o = orient_segment(s1, sp, orientation);
	o2 = vv_oriented_side(sv, v2->site(), s1, s1o, sp);
	o3 = vv_oriented_side(sv, s1, v3->site(), s1o, sp);
	CGAL_assertion( o2 != ON_ORIENTED_BOUNDARY );
	CGAL_assertion( o3 != ON_ORIENTED_BOUNDARY );
	if ( o2 == ON_POSITIVE_SIDE && o3 == ON_NEGATIVE_SIDE ) {
	  return *ec;
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

#endif // CGAL_VORONOI_DIAGRAM_2_SEGMENT_DELAUNAY_GRAPH_NEAREST_SITE_2_H
