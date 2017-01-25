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

#ifndef CGAL_VORONOI_DIAGRAM_2_HALFEDGE_2_H
#define CGAL_VORONOI_DIAGRAM_2_HALFEDGE_2_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Finder_classes.h>
#include <CGAL/Triangulation_utils_2.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

template<class VDA>
class Halfedge
{
 private:
  typedef Halfedge<VDA>                          Self;
  typedef Triangulation_cw_ccw_2                 CW_CCW_2;

  typedef typename VDA::Delaunay_graph           DG;
  typedef typename DG::Face_handle               Delaunay_face_handle;
  typedef typename DG::Vertex_circulator         Delaunay_vertex_circulator;
  typedef typename DG::Edge_circulator           Delaunay_edge_circulator;
  typedef typename VDA::Accessor::Edge_rejector  Edge_rejector;
  typedef typename VDA::Accessor::Face_rejector  Face_rejector;

 private:
  Self find_next_1D() const {
    Delaunay_vertex_circulator vc = vda_->dual().incident_vertices(v1_);

    while ( v2_ == vc ) { ++vc; }

    if ( vda_->dual().is_infinite(vc) ) { return *this; }

    return Self(vda_, v1_, vc);
  }

  void find_next(const Delaunay_face_handle& f, int i,
		 Delaunay_face_handle& fnext, int& inext) const {
    Find_next_halfedge<VDA>()(vda_, f, i, fnext, inext);
  }

  void find_opposite(const Delaunay_face_handle& f, int i,
		     Delaunay_face_handle& fopp, int& iopp) const {
    Find_opposite_halfedge<VDA>()(vda_, f, i, fopp, iopp);
  }

public:
  typedef typename VDA::Vertex                    Vertex;
  typedef typename VDA::Face                      Face;
  typedef typename VDA::Vertex_handle             Vertex_handle;
  typedef typename VDA::Face_handle               Face_handle;
  typedef typename VDA::Halfedge_handle           Halfedge_handle;
  typedef typename VDA::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
  typedef typename VDA::Delaunay_graph            Delaunay_graph;
  typedef typename Delaunay_graph::Edge           Delaunay_edge;
  typedef typename Delaunay_graph::Vertex_handle  Delaunay_vertex_handle;

  // CONSTRUCTORS
  //-------------
  Halfedge(const VDA* vda = NULL)
    : vda_(vda), f_(Delaunay_face_handle()), i_(-1),
      v1_(Delaunay_vertex_handle()), v2_(Delaunay_vertex_handle()) {}
  Halfedge(const VDA* vda, Delaunay_face_handle f, int i)
    : vda_(vda), f_(f), i_(i), v1_(Delaunay_vertex_handle()),
      v2_(Delaunay_vertex_handle())
  {
    CGAL_precondition( vda_->dual().dimension() == 2 );
    CGAL_precondition( !vda_->edge_rejector()(vda_->dual(), f_, i_) );
  }

  Halfedge(const VDA* vda, Delaunay_vertex_handle v1,
	   Delaunay_vertex_handle v2)
    : vda_(vda), f_(Delaunay_face_handle()), i_(-2), v1_(v1), v2_(v2)
  {
    CGAL_precondition( vda_->dual().dimension() == 1 );
    CGAL_precondition( !vda_->dual().is_infinite(v1_) );
  }

  // ACCESS TO NEIGHBORING HALFEDGES
  //--------------------------------
  Halfedge_handle opposite() const {
    if ( vda_->dual().dimension() == 1 ) {
      return Halfedge_handle( Self(vda_, v2_, v1_) );
    }

    int cw_i = CW_CCW_2::cw(i_);
    if (  vda_->face_rejector()(vda_->dual(), f_->vertex(cw_i))  ) {
      Delaunay_face_handle fopp;
      int iopp;
      find_opposite(f_, i_, fopp, iopp); //equivalent to: twin().next().twin();
      return Halfedge_handle( Self(vda_, fopp, iopp) );
    } else {
      int i_mirror = vda_->dual().tds().mirror_index(f_, i_);
      return
	Halfedge_handle( Self(vda_, f_->neighbor(i_), i_mirror) );
    }
  }

  inline Halfedge_handle twin() const {
    return opposite();
  }

  Halfedge_handle next() const {
    if ( vda_->dual().dimension() == 1 ) {
      return Halfedge_handle( find_next_1D() );
    }

    // if I want to return all edges and not just the finite ones,
    // replace the do-while loop by the following statements:
    //          find_next(f_, i_, f, i);
    Delaunay_face_handle f = f_, fnext;
    int i = i_, inext;
    do {
      find_next(f, i, fnext, inext);
      f = fnext;
      i = inext;
    } while ( vda_->dual().is_infinite(f, i) );
    return Halfedge_handle( Self(vda_, f, i) );
  }

  Halfedge_handle previous() const
  {
    if ( vda_->dual().dimension() == 1 ) {
      return Halfedge_handle( find_next_1D() );
    }

    Delaunay_face_handle f, fprev = f_;
    int iprev = i_, i;
    
    // if I want to return also infinite edges replace the test in
    // the while loop by the following test (i.e., should omit the
    // testing for infinity):
    //           vda_->edge_rejector()(vda_->dual(), f, i)
    do {
      f = fprev->neighbor(iprev);
      int i_mirror = vda_->dual().tds().mirror_index(fprev, iprev);
      i = CW_CCW_2::ccw( i_mirror );
      fprev = f;
      iprev = i;
    } while ( vda_->edge_rejector()(vda_->dual(), f, i) ||
	      vda_->dual().is_infinite(f, i) );

    return Halfedge_handle( Self(vda_, f, i) );
  }

  Vertex_handle source() const {
    CGAL_precondition( has_source() );
    return opposite()->target();
  }

  Vertex_handle target() const {
    CGAL_precondition( has_target() );

    Delaunay_face_handle fvalid = Find_valid_vertex<VDA>()(vda_, f_);
    CGAL_assertion( !vda_->dual().is_infinite(fvalid) );
    return Vertex_handle( Vertex(vda_, fvalid) );
  }

  Face_handle face() const {
    if ( vda_->dual().dimension() == 1 ) {
      return Face_handle( Face(vda_, v1_) );
    }
    Face f(vda_, f_->vertex( CW_CCW_2::ccw(i_) ));
    return Face_handle(f);
    //      Face_handle(   Face(vda_, f_->vertex( CW_CCW_2::ccw(i_) ))   );
  }

  Ccb_halfedge_circulator ccb() const {
    return Ccb_halfedge_circulator( *this );
  }

  // PREDICATES
  //-----------
  bool has_source() const {
    return opposite()->has_target();
  }

  bool has_target() const {
    if ( vda_->dual().dimension() == 1 ) { return false; }
    return !vda_->dual().is_infinite(f_);
  }

  bool is_unbounded() const {
    return !has_source() || !has_target();
  }

  bool is_bisector() const {
    return !has_source() && !has_target();
  }

  bool is_segment() const {
    return has_source() && has_target();
  }

  bool is_ray() const {
    return
      ( has_source() && !has_target() ) ||
      ( !has_source() && has_target() );
  }

  // ACCESS TO DUAL DEFINING VERTICES
  //---------------------------------
  Delaunay_vertex_handle up() const {
    if ( vda_->dual().dimension() == 1 ) { return v1_; }
    return f_->vertex( CW_CCW_2::ccw(i_) );
  }

  Delaunay_vertex_handle down() const {
    if ( vda_->dual().dimension() == 1 ) { return v2_; }
    return f_->vertex( CW_CCW_2::cw(i_) );
  }

  Delaunay_vertex_handle left() const {
    CGAL_precondition( has_source() );
    return vda_->dual().tds().mirror_vertex(f_,i_);
  }

  Delaunay_vertex_handle right() const {
    CGAL_precondition( has_target() );
    return f_->vertex( i_ );
  }


  // DUAL FEATURE
  //-------------
  Delaunay_edge dual() const {
    if ( vda_->dual().dimension() == 1 ) {
      Delaunay_edge_circulator ec;
      if ( vda_->dual().is_infinite(v1_) ) {
	CGAL_assertion( !vda_->dual().is_infinite(v2_) );
	ec = vda_->dual().incident_edges(v2_);
      } else {
	ec = vda_->dual().incident_edges(v1_);
      }
      do {
	Delaunay_edge e = *ec;
	if ( (e.first->vertex(0) == v1_ && e.first->vertex(1) == v2_) ||
	     (e.first->vertex(0) == v2_ && e.first->vertex(1) == v1_) ) {
	  return e;
	}
	++ec;
      } while ( true );
    } else {
      return Delaunay_edge(f_, i_);
    }
  }

  // VALIDITY TESTING
  //-----------------
  bool is_valid() const {
    if ( vda_ == NULL ) { return true; }

    bool valid = true;

    if ( vda_->dual().dimension() == 1 ) {
      valid = valid && v1_ != Delaunay_vertex_handle();
      valid = valid && v2_ != Delaunay_vertex_handle();
    } else {
      valid = valid && !vda_->edge_rejector()(vda_->dual(), f_, i_);

      Delaunay_vertex_handle v = f_->vertex( CW_CCW_2::ccw(i_) );

      valid = valid && !vda_->face_rejector()(vda_->dual(), v);
    }

    Halfedge_handle h_this(*this);

    valid = valid && opposite()->opposite() == h_this;

    if ( has_source() ) {
      valid = valid && source()->is_valid();
      valid = valid && source()->is_incident_edge( h_this );
    }

    if ( has_target() ) {
      valid = valid && target()->is_valid();
      valid = valid && target()->is_incident_edge( h_this );
    }

    valid = valid && next()->previous() == h_this;
    valid = valid && previous()->next() == h_this;
    return valid;
  }

  // COMPARISON OPERATORS
  //---------------------
  bool operator==(const Self& other) const {
    if ( vda_ == NULL ) { return other.vda_ == NULL; }
    if ( other.vda_ == NULL ) { return vda_ == NULL; }

    if ( vda_->dual().dimension() == 1 ) {
      return ( vda_ == other.vda_ && v1_ == other.v1_ && v2_ == other.v2_ );
    } else {
      return ( vda_ == other.vda_ && f_ == other.f_ && i_ == other.i_ );
    }
  }
  
  bool operator!=(const Self& other) const {
    return !((*this) == other);
  }

  bool operator<(const Self& other) const {
    if ( vda_ == NULL ) { return other.vda_ != NULL; }
    if ( other.vda_ == NULL ) { return false; }

    if ( vda_ != other.vda_ ) { return vda_ < other.vda_; }
    if ( vda_->dual().dimension() == 1 ) {
      if ( v1_ != other.v1_ ) { return v1_ < other.v1_; }
      return v2_ < other.v2_;
    } else {
      if ( f_ != other.f_ ) { return f_ < other.f_; }
      return i_ < other.i_;
    }
  }

 private:
  const VDA* vda_;
  Delaunay_face_handle f_;
  int i_;
  Delaunay_vertex_handle v1_, v2_;
};


} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

//=========================================================================
//=========================================================================

#endif // CGAL_VORONOI_DIAGRAM_2_HALFEDGE_2_H
