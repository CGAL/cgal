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

#ifndef CGAL_VORONOI_DIAGRAM_2_FACE_H
#define CGAL_VORONOI_DIAGRAM_2_FACE_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Accessor.h>
#include <CGAL/Triangulation_utils_2.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

template<class VDA>
class Face
{
 private:
  typedef typename Accessor<VDA>::Non_degenerate_faces_iterator
  Non_degenerate_faces_iterator;

 private:
  typedef Face<VDA>                              Self;
  typedef typename VDA::Delaunay_graph           DG;
  typedef typename DG::Edge_circulator           Delaunay_edge_circulator;
  typedef typename DG::Vertex_circulator         Delaunay_vertex_circulator;

  typedef Triangulation_cw_ccw_2                 CW_CCW_2;

 public:
  typedef typename VDA::Delaunay_graph           Delaunay_graph;
  typedef typename VDA::Halfedge                 Halfedge;
  typedef typename VDA::Vertex                   Vertex;
  typedef typename VDA::Halfedge_handle          Halfedge_handle;
  typedef typename VDA::Vertex_handle            Vertex_handle;
  typedef typename VDA::Face_handle              Face_handle;
  typedef typename VDA::Ccb_halfedge_circulator  Ccb_halfedge_circulator;
  typedef typename Delaunay_graph::Vertex_handle Delaunay_vertex_handle;

 public:

  // CONSTRUCTORS
  //-------------
  Face(const VDA* vda = NULL) : vda_(vda) {}
  Face(const VDA* vda, Delaunay_vertex_handle v) : vda_(vda), v_(v)
  {
    //    CGAL_precondition( !vda_->face_rejector()(v_) );
  }

  Face(const VDA* vda, const Non_degenerate_faces_iterator& it)
    : vda_(vda), v_(it.base())
  {
    //    CGAL_precondition( !vda_->face_rejector()(v_) );
  }

  // ACCESS TO NEIGHBORING OBJECTS
  //------------------------------
  Halfedge_handle halfedge() const
  {
    CGAL_precondition( vda_->dual().dimension() > 0 );
    if ( vda_->dual().dimension() == 1 ) {
      Delaunay_vertex_circulator vc;
      vc = vda_->dual().incident_vertices(v_);
      while ( vda_->dual().is_infinite(vc) ) { ++vc; }
      Delaunay_vertex_handle vv(vc);
      return Halfedge_handle( Halfedge(vda_, v_, vv) );
    }

    // the edge circulator gives edges that have v_ as their target
    Delaunay_edge_circulator ec = vda_->dual().incident_edges(v_);
    CGAL_assertion_code( Delaunay_edge_circulator ec_start = ec );

    // if I want to return also infinite edges replace the test in
    // the while loop by the following test (i.e., should omit the
    // testing for infinity):
    //           vda_->edge_rejector()(vda_->dual(), ec)
    while ( vda_->edge_rejector()(vda_->dual(), ec) ||
	    vda_->dual().is_infinite(ec) ) {
      ++ec;
      CGAL_assertion( ec != ec_start );
    }
    CGAL_assertion(ec->first->vertex( CW_CCW_2::cw(ec->second) ) == v_);

    int i_mirror = vda_->dual().tds().mirror_index(ec->first, ec->second);

#if !defined(CGAL_NO_ASSERTIONS) && !defined(NDEBUG)
    Halfedge h(vda_, ec->first->neighbor(ec->second), i_mirror);
    Face_handle f_this(*this);
    CGAL_assertion( h.face() == f_this );
#endif

    return
      Halfedge_handle( 
		      Halfedge(vda_, ec->first->neighbor(ec->second), i_mirror)
		      );
  }

  Ccb_halfedge_circulator ccb() const {
    return Ccb_halfedge_circulator( *halfedge() );
  }

  Ccb_halfedge_circulator outer_ccb() const { return ccb(); }
  // PREDICATES
  //-----------
  bool is_unbounded() const {
    if ( vda_->dual().dimension() < 2 ) { return true; }

    Delaunay_vertex_circulator vc = vda_->dual().incident_vertices(v_);
    Delaunay_vertex_circulator vc_start = vc;
    do {
      if ( vda_->dual().is_infinite(vc) ) { return true; }
      ++vc;
    } while ( vc != vc_start );
    return false;
  }

  bool is_halfedge_on_ccb(const Halfedge_handle& he) const {
    Ccb_halfedge_circulator hc_start = ccb();
    Ccb_halfedge_circulator hc = hc_start;
    do {
      if ( he == *hc ) { return true; }
      hc++;
    } while ( hc != hc_start );
    return false;
  }

  // DUAL FEATURE
  //-------------
  const Delaunay_vertex_handle& dual() const { return v_; }

  // VALIDITY TESTING
  //-----------------
  bool is_valid() const {
    if ( vda_ == NULL ) { return true; }

    if ( vda_->dual().dimension() < 1 ) { return true; }

    bool valid = !vda_->face_rejector()(vda_->dual(), v_);

    valid = valid && !vda_->edge_rejector()(vda_->dual(), halfedge()->dual());

    Ccb_halfedge_circulator hc = ccb();
    Ccb_halfedge_circulator hc_start = hc;
    Face_handle f_this(*this);
    do {
      valid = valid && hc->face() == f_this;
      valid = valid && !vda_->edge_rejector()( vda_->dual(), hc->dual() );
      hc++;
    } while ( hc != hc_start );

    return valid;
  }

  // COMPARISON OPERATORS
  //---------------------
  bool operator==(const Self& other) const {
    if ( vda_ == NULL ) { return other.vda_ == NULL; }
    if ( other.vda_ == NULL ) { return vda_ == NULL; }
    return ( vda_ == other.vda_ && v_ == other.v_ );
  }

  bool operator!=(const Self& other) const {
    return !((*this) == other);
  }

  bool operator<(const Self& other) const {
    if ( vda_ == NULL ) { return other.vda_ != NULL; }
    if ( other.vda_ == NULL ) { return false; }
    if ( vda_ != other.vda_ ) { return vda_ < other.vda_; }
    return v_ < other.v_;
  }

private:
  const VDA* vda_;
  Delaunay_vertex_handle v_;
};

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL


#endif // CGAL_VORONOI_DIAGRAM_2_FACE_H
