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

#ifndef CGAL_VORONOI_DIAGRAM_2_DELAUNAY_TRIANGULATION_NEAREST_SITE_2_H
#define CGAL_VORONOI_DIAGRAM_2_DELAUNAY_TRIANGULATION_NEAREST_SITE_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Triangulation_utils_2.h>

#include <boost/variant.hpp>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================


template<class DG>
class Delaunay_triangulation_nearest_site_2
{
 public:
  typedef DG                                          Delaunay_graph;
  typedef typename DG::Geom_traits::Point_2           Point_2;

 private:
  typedef Triangulation_cw_ccw_2                      CW_CCW_2;

  typedef typename Delaunay_graph::Geom_traits        Geom_traits;
  typedef typename Delaunay_graph::Vertex_handle      Vertex_handle;
  typedef typename Delaunay_graph::Face_handle        Face_handle;
  typedef typename Delaunay_graph::Edge               Edge;
  typedef typename Delaunay_graph::Vertex_circulator  Vertex_circulator;
  typedef typename Delaunay_graph::Face_circulator    Face_circulator;
  typedef typename Delaunay_graph::Edge_circulator    Edge_circulator;

 public:
  typedef boost::variant<Vertex_handle,Edge,Face_handle>  result_type;

  result_type operator()(const Delaunay_graph& dg, const Point_2& p) const {
    CGAL_precondition( dg.dimension() >= 0 );

    typename Geom_traits::Compare_distance_2 compare_distance =
      dg.geom_traits().compare_distance_2_object();

    Vertex_handle v = dg.nearest_vertex(p);

    if ( dg.dimension() == 0 ) {
      return v;
    }

    if ( dg.dimension() == 1 ) {
      Edge_circulator ec = dg.incident_edges(v);
      Edge_circulator ec_start = ec;
      Comparison_result cr;

      do {
	Edge e = *ec;
	Vertex_handle v1 = e.first->vertex(CW_CCW_2::ccw(e.second));
	Vertex_handle v2 = e.first->vertex(CW_CCW_2::cw(e.second) );

	if ( v == v1 ) {
	  if ( !dg.is_infinite(v2) ) {
	    cr = compare_distance(p, v2->point(), v->point());
	    CGAL_assertion( cr != SMALLER );
	    if ( cr == EQUAL ) {
	      return e;
	    }
	  }
	} else {
	  CGAL_assertion( v == v2 );
	  if ( !dg.is_infinite(v1) ) {
	    cr = compare_distance(p, v1->point(), v->point());
	    CGAL_assertion( cr != SMALLER );
	    if ( cr == EQUAL ) {
	      return e;
	    }
	  }
	}
	++ec;
      } while ( ec != ec_start );

      return v;
    }

    CGAL_assertion( dg.dimension() == 2 );

    Face_circulator fc_start = dg.incident_faces(v);
    Face_circulator fc = fc_start;

    // first check if the point lies on a Voronoi vertex
    do {
      int index = fc->index(v);
      Vertex_handle v1 = fc->vertex( CW_CCW_2::ccw(index) );
      Vertex_handle v2 = fc->vertex( CW_CCW_2::cw(index)  );

      Comparison_result cr1 = LARGER, cr2 = LARGER;

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	cr1 = compare_distance(p, v1->point(), v->point());
      }
      if ( !dg.is_infinite(v2) ) {
	cr2 = compare_distance(p, v2->point(), v->point());
      }

      CGAL_assertion( cr1 != SMALLER );
      CGAL_assertion( cr2 != SMALLER );

      if ( cr1 == EQUAL && cr2 == EQUAL ) {
	return Face_handle(fc);
      }

      ++fc;
    } while ( fc != fc_start );

    // now check if the point lies on a Voronoi edge
    fc_start = dg.incident_faces(v);
    fc = fc_start;
    do {
      int index = fc->index(v);
      Vertex_handle v1 = fc->vertex( CW_CCW_2::ccw(index) );
      Vertex_handle v2 = fc->vertex( CW_CCW_2::cw(index)  );

      Comparison_result cr1 = LARGER, cr2 = LARGER;

      // do the generic check now
      if ( !dg.is_infinite(v1) ) {
	cr1 = compare_distance(p, v1->point(), v->point());
      }
      if ( !dg.is_infinite(v2) ) {
	cr2 = compare_distance(p, v2->point(), v->point());
      }

      CGAL_assertion( cr1 != SMALLER );
      CGAL_assertion( cr2 != SMALLER );
      CGAL_assertion( cr1 != EQUAL || cr2 != EQUAL );

      if ( cr1 == EQUAL ) {
	Face_handle f(fc);
	return Edge(f, CW_CCW_2::cw(index) );
      } else if ( cr2 == EQUAL ) {
	Face_handle f(fc);
	return Edge(f, CW_CCW_2::ccw(index) );
      }

      ++fc;
    } while ( fc != fc_start );

    return v;
  }
};

//=========================================================================
//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_DELAUNAY_TRIANGULATION_NEAREST_SITE_2_H
