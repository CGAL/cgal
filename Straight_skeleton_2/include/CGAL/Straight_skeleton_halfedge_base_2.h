// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H
#define CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/tags.h>
#include <CGAL/enum.h>

namespace CGAL {

template < class Refs >
class Straight_skeleton_halfedge_base_base_2
{
public:

  typedef Refs     HalfedgeDS;
  typedef Tag_true Supports_halfedge_prev;
  typedef Tag_true Supports_halfedge_vertex;
  typedef Tag_true Supports_halfedge_face;

  typedef typename Refs::Vertex_handle         Vertex_handle;
  typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
  typedef typename Refs::Halfedge_handle       Halfedge_handle;
  typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Refs::Face_handle           Face_handle;
  typedef typename Refs::Face_const_handle     Face_const_handle;
  typedef typename Refs::Vertex                Vertex;
  typedef typename Refs::Face                  Face;

  typedef Straight_skeleton_halfedge_base_base_2<Refs> Base_base ;

protected:

  Straight_skeleton_halfedge_base_base_2() : mF(Face_handle()), mID(-1), mSlope(ZERO) {}

  Straight_skeleton_halfedge_base_base_2( int aID ) : mF(Face_handle()), mID(aID), mSlope(ZERO) {}

  Straight_skeleton_halfedge_base_base_2( int aID, Sign aSlope ) : mF(Face_handle()), mID(aID), mSlope(aSlope) {}

public:

  int id() const { return mID ; }

  bool is_bisector() const
  {
    return !this->is_border() && !this->opposite()->is_border() ;
  }

  bool is_inner_bisector() const
  {
    return !this->vertex()->is_contour() && !this->opposite()->vertex()->is_contour();
  }

  bool has_null_segment() const { return this->vertex()->has_null_point() ; }

  bool has_infinite_time() const { return this->vertex()->has_infinite_time() ; }

  Halfedge_const_handle defining_contour_edge() const { return this->face()->halfedge() ; }
  Halfedge_handle       defining_contour_edge()       { return this->face()->halfedge() ; }

  Halfedge_handle       opposite()       { return mOpp;}
  Halfedge_const_handle opposite() const { return mOpp;}
  Halfedge_handle       next    ()       { return mNxt;}
  Halfedge_const_handle next    () const { return mNxt;}
  Halfedge_handle       prev    ()       { return mPrv; }
  Halfedge_const_handle prev    () const { return mPrv; }
  Vertex_handle         vertex  ()       { return mV; }
  Vertex_const_handle   vertex  () const { return mV; }
  Face_handle           face    ()       { return mF; }
  Face_const_handle     face    () const { return mF; }

  Sign slope() const { return mSlope ; }

  bool is_border() const { return mF == Face_handle();}

  void set_opposite( Halfedge_handle h) { mOpp = h; }
  void set_next    ( Halfedge_handle h) { mNxt = h; }
  void set_prev    ( Halfedge_handle h) { mPrv = h; }
  void set_vertex  ( Vertex_handle   w) { mV   = w; }
  void set_face    ( Face_handle     g) { mF   = g; }

  void set_slope( Sign aSlope ) { mSlope = aSlope ; }

  void reset_id ( int aID ) { mID = aID ; }

private:

  Halfedge_handle  mOpp;
  Halfedge_handle  mNxt;
  Halfedge_handle  mPrv;
  Vertex_handle    mV;
  Face_handle      mF;
  int              mID ;
  Sign             mSlope ;
};

template < class Refs >
class Straight_skeleton_halfedge_base_2 : public Straight_skeleton_halfedge_base_base_2<Refs>
{
public:

  typedef typename Refs::Vertex_handle   Vertex_handle;
  typedef typename Refs::Halfedge_handle Halfedge_handle;
  typedef typename Refs::Face_handle     Face_handle;

  typedef Straight_skeleton_halfedge_base_base_2<Refs> Base_base ;
  typedef Straight_skeleton_halfedge_base_2<Refs>      Base ;

  Straight_skeleton_halfedge_base_2() {}

  Straight_skeleton_halfedge_base_2( int aID ) : Base_base(aID) {}

  Straight_skeleton_halfedge_base_2( int aID, Sign aSlope ) : Base_base(aID,aSlope) {}

  void set_opposite( Halfedge_handle h )  { Base_base::opposite(h)  ; }
  void set_next    ( Halfedge_handle h )  { Base_base::set_next(h)  ; }
  void set_prev    ( Halfedge_handle h )  { Base_base::set_prev(h)  ; }
  void set_vertex  ( Vertex_handle   w )  { Base_base::set_vertex(w); }
  void set_face    ( Face_handle     g )  { Base_base::set_face(g)  ; }
  void set_slope   ( Sign            s )  { Base_base::set_slope(s) ; }
  void reset_id    ( int             i )  { Base_base::reset_id(i) ; }

} ;

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H //
// EOF //

