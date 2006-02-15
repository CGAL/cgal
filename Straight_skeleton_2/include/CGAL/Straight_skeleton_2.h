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
  typedef typename Base::Halfedge_iterator Halfedge_iterator  ;

  Straight_skeleton_2() {}

/*
private:

  Face_handle     faces_push_back    ( Face const& aF ) ;
  Vertex_handle   vertices_push_back ( Vertex const& aV ) ;
  Halfedge_handle edges_push_back    ( Halfedge const& aA, Halfedge const& aB )  ;

  void vertices_erase ( Vertex_iterator   first, Vertex_iterator   last )  ;
  void edges_erase    ( Halfedge_handle h )  ;
*/
};

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_2_H //
// EOF //

