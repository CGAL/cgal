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
// file          : include/CGAL/Straight_skeleton_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================

#ifndef CGAL_STRAIGHT_SKELETON_2_H
#define CGAL_STRAIGHT_SKELETON_2_H 1

#ifndef CGAL_STRAIGHT_SKELETON_ITEMS_2_H
#include <CGAL/Straight_skeleton_items_2.h>
#endif

#ifndef CGAL_HALFEDGEDS_DEFAULT_H
#include <CGAL/HalfedgeDS_default.h>
#endif

CGAL_BEGIN_NAMESPACE

template<  class Traits_
         , class Items_ = Straight_skeleton_items_2
         , class Alloc_ = CGAL_ALLOCATOR(int)
        >
class Straight_skeleton_2 : public CGAL_HALFEDGEDS_DEFAULT <Traits_,Items_,Alloc_>
{
public :

  typedef Traits_ Rep ;
  
  typedef CGAL_HALFEDGEDS_DEFAULT <Traits_,Items_,Alloc_> Base ;
  
  typedef typename Base::Vertex_base     Vertex ;
  typedef typename Base::Halfedge_base   Halfedge ;
  typedef typename Base::Face_base       Face ;
  typedef typename Base::Vertex_handle   Vertex_handle ;
  typedef typename Base::Halfedge_handle Halfedge_handle  ;
  typedef typename Base::Face_handle     Face_handle  ;
  typedef typename Base::Vertex_iterator Vertex_iterator ;
  
  Straight_skeleton_2() {}
  
private:
 
  Face_handle     faces_push_back    ( Face const& aF ) ;
  Vertex_handle   vertices_push_back ( Vertex const& aV ) ;
  Halfedge_handle edges_push_back    ( Halfedge const& aA, Halfedge const& aB )  ;
  
  void vertices_erase ( Vertex_iterator   first, Vertex_iterator   last )  ;
  void edges_erase    ( Halfedge_iterator first, Halfedge_iterator last )  ;
    
  /*
  Vertex_iterator   vertices_begin () ;
  Vertex_iterator   vertices_end   () ;
  Halfedge_iterator halfedges_begin() ;
  Halfedge_iterator halfedges_end  () ;
  */
  
};

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_2_H //
// EOF //
 
