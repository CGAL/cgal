// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Straight_skeleton_halfedge_base_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================

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
    return !is_border() && !opposite()->is_border() ;
  }
  
  bool is_contour_bisector() const
  {
    return vertex()->is_border() || opposite()->vertex()->is_border();
  }
  
  Halfedge_const_handle defining_border() const { return face()->halfedge() ; }
  Halfedge_handle       defining_border()       { return face()->halfedge() ; }
  
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
 
