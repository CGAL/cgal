// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H
#define CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H 1

#ifndef CGAL_HALFEDGEDS_HALFEDGE_BASE_H
#include <CGAL/HalfedgeDS_halfedge_base.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class Refs, class S >
class Straight_skeleton_halfedge_base_base_2
  : public HalfedgeDS_halfedge_base<Refs, Tag_true, Tag_true, Tag_true >
{
public:

  typedef HalfedgeDS_halfedge_base<Refs, Tag_true, Tag_true, Tag_true> HDSBase ;

  typedef typename HDSBase::Base HBase ;

  typedef typename HDSBase::Base_base HBase_base ;

  typedef S Segment_2;

  typedef typename Refs::Halfedge_handle       Halfedge_handle ;
  typedef typename Refs::Halfedge_const_handle Halfedge_const_handle ;
public:

  Straight_skeleton_halfedge_base_base_2() : mID(-1) {}
  Straight_skeleton_halfedge_base_base_2 ( int aID ) : mID(aID) {}

  int id() const { return mID ; }

  bool is_bisector() const
  {
    return !HBase::is_border() && !HBase::opposite()->is_border() ;
  }

  bool is_inner_bisector() const
  {
    return !vertex()->is_contour() && !opposite()->vertex()->is_contour();
  }

  Halfedge_const_handle defining_contour_edge() const { return HBase::face()->halfedge() ; }
  Halfedge_handle       defining_contour_edge()       { return HBase::face()->halfedge() ; }

  void  set_opposite( Halfedge_handle h)  { HBase::set_opposite(h);}


private:
  int mID ;
};

template < class Refs, class S >
class Straight_skeleton_halfedge_base_2
  : public Straight_skeleton_halfedge_base_base_2 < Refs, S >
{
public:

  typedef Straight_skeleton_halfedge_base_base_2 < Refs, S > SSBase ;

  typedef typename SSBase::HBase_base HBase_base ;

  typedef typename SSBase::Base      HBase ;

  typedef typename SSBase::Segment_2 Segment_2;
  typedef typename HBase::Halfedge_handle Halfedge_handle ;
  typedef typename HBase::Vertex_handle   Vertex_handle ;
  typedef typename HBase::Face_handle     Face_handle ;

public:

  Straight_skeleton_halfedge_base_2() {}
  Straight_skeleton_halfedge_base_2( int aID ) : SSBase(aID) {}

  void  set_opposite( Halfedge_handle h)  { HBase_base::set_opposite(h);}

protected:

  void set_prev  ( Halfedge_handle h ) { HBase_base::set_prev(h) ; }
  void set_vertex( Vertex_handle   w ) { HBase_base::set_vertex(w); }
  void set_face  ( Face_handle     g ) { HBase_base::set_face(g) ; }

};
CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H //
// EOF //

