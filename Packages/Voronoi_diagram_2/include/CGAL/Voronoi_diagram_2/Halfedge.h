// Copyright (c) 2005 University of Crete (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_HALFEDGE_2_H
#define CGAL_VORONOI_DIAGRAM_2_HALFEDGE_2_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Finder_classes.h>
#include <CGAL/Triangulation_utils_2.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

template<class VDA>
class Halfedge
{
 private:
  typedef Halfedge<VDA>                         Self;
  typedef Triangulation_cw_ccw_2                CW_CCW_2;

  typedef typename VDA::Dual_face_handle        Dual_face_handle;
  typedef typename VDA::Dual_vertex_handle      Dual_vertex_handle;
  typedef typename VDA::Edge_degeneracy_tester  Edge_degeneracy_tester;
  typedef typename VDA::Face_degeneracy_tester  Face_degeneracy_tester;

 private:
  void find_next(const Dual_face_handle& f, int i,
		 Dual_face_handle& fnext, int& inext) const
  {
    Find_next_halfedge<VDA>()(vda_, f, i, fnext, inext);
  }

  void find_opposite(const Dual_face_handle& f, int i,
		     Dual_face_handle& fopp, int& iopp) const 
  {
    Find_opposite_halfedge<VDA>()(vda_, f, i, fopp, iopp);
  }

public:
  typedef typename VDA::Vertex                   Vertex;
  typedef typename VDA::Face                     Face;
  typedef typename VDA::Vertex_handle            Vertex_handle;
  typedef typename VDA::Face_handle              Face_handle;
  typedef typename VDA::Halfedge_handle          Halfedge_handle;
  typedef typename VDA::Ccb_halfedge_circulator  Ccb_halfedge_circulator;


  Halfedge(const VDA* vda = NULL)
    : vda_(vda), f_(Dual_face_handle()), i_(-1) {}
  Halfedge(const VDA* vda, Dual_face_handle f, int i)
    : vda_(vda), f_(f), i_(i)
  {
#ifndef CGAL_NO_PRECONDITIONS
    CGAL_precondition( !vda_->edge_tester()(vda_->dual(), f_, i_) );
    Dual_vertex_handle v = f_->vertex( CW_CCW_2::ccw(i_) );
#if 0
    // the following test prohibit creating a halfedge which is
    // temporarily invalid...
    CGAL_precondition( !vda_->face_tester()(v) );
#endif
#endif
  }

  Halfedge_handle opposite() const {
    int cw_i = CW_CCW_2::cw(i_);
    if ( vda_->face_tester()(vda_->dual(), vda_->edge_tester(),
			     f_->vertex(cw_i)) ) {
      Dual_face_handle fopp;
      int iopp;
      find_opposite(f_, i_, fopp, iopp); //equivalent to: twin().next().twin();
      return Halfedge_handle( Self(vda_, fopp, iopp) );
    } else {
      int i_mirror = vda_->dual().tds().mirror_index(f_, i_);
      return
	Halfedge_handle( Self(vda_, f_->neighbor(i_), i_mirror) );
    }
  }

  Halfedge_handle twin() const {
    return opposite();
  }

  Halfedge_handle next_halfedge() const {
    return next();
  }

  Halfedge_handle next() const {
    // if I want to return all edges and not just the finite ones,
    // replace the do-while loop by the following statements:
    //          find_next(f_, i_, f, i);
    Dual_face_handle f = f_, fnext;
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
    Dual_face_handle f, fprev = f_;
    int iprev = i_, i;
    
    // if I want to return also infinite edges replace the test in
    // the while loop by the following test (i.e., should omit the
    // testing for infinity):
    //           vda_->edge_tester()(vda_->dual(), f, i)
    do {
      f = fprev->neighbor(iprev);
      int i_mirror = vda_->dual().tds().mirror_index(fprev, iprev);
      i = CW_CCW_2::ccw( i_mirror );
      fprev = f;
      iprev = i;
    } while ( vda_->edge_tester()(vda_->dual(), f, i) ||
	      vda_->dual().is_infinite(f, i) );

    return Halfedge_handle( Self(vda_, f, i) );
  }

  Ccb_halfedge_circulator ccb() const {
    return Ccb_halfedge_circulator( Halfedge_handle(*this) );
  }

  bool operator==(const Self& other) const {
    if ( vda_ == NULL ) { return other.vda_ == NULL; }
    if ( other.vda_ == NULL ) { return vda_ == NULL; }
    return ( vda_ == other.vda_ && f_ == other.f_ && i_ == other.i_ );
  }
  
  bool operator!=(const Self& other) const {
    return !((*this) == other);
  }

  Vertex_handle source() const {
    return opposite()->vertex();
  }

  Vertex_handle target() const {
    return vertex();
  }

  Vertex_handle vertex() const {
    Dual_face_handle fvalid = Find_valid_vertex<VDA>()(vda_, f_);
    CGAL_assertion( !vda_->dual().is_infinite(fvalid) );
    return Vertex_handle( Vertex(vda_, fvalid) );
  }

  Face_handle face() const {
    return
      Face_handle(   Face(vda_, f_->vertex( CW_CCW_2::ccw(i_) ))   );
  }

  // this is temporary or maybe not
  typename VDA::Dual_edge dual_edge() const {
    return typename VDA::Dual_edge(f_, i_);
  }

  bool is_valid() const {
    if ( vda_ == NULL ) { return true; }

    bool valid = !vda_->edge_tester()(vda_->dual(), f_, i_);

    Dual_vertex_handle v = f_->vertex( CW_CCW_2::ccw(i_) );

    valid = valid && !vda_->face_tester()(vda_->dual(),vda_->edge_tester(),v);

    Halfedge_handle h_this(*this);

    valid = valid && opposite()->opposite() == h_this;

    valid = valid && source()->is_valid();
    valid = valid && target()->is_valid();

    valid = valid && source()->is_incident_edge( h_this );
    valid = valid && target()->is_incident_edge( h_this );

    valid = valid && next()->previous() == h_this;
    valid = valid && previous()->next() == h_this;
    return valid;
  }

private:
  const VDA* vda_;
  Dual_face_handle f_;
  int i_;
};

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_HALFEDGE_2_H
