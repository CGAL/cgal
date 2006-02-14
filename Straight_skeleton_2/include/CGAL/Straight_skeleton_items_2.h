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
// file          : include/CGAL/Straight_skeleton_items_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================

#ifndef CGAL_STRAIGHT_SKELETON_ITEMS_2_H
#define CGAL_STRAIGHT_SKELETON_ITEMS_2_H 1

#ifndef CGAL_STRAIGHT_SKELETON_VERTEX_BASE_2_H
#include <CGAL/Straight_skeleton_vertex_base_2.h>
#endif

#ifndef CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H
#include <CGAL/Straight_skeleton_halfedge_base_2.h>
#endif

#ifndef CGAL_HALFEDGEDS_FACE_BASE_H
#include <CGAL/HalfedgeDS_face_base.h>
#endif

CGAL_BEGIN_NAMESPACE

class Straight_skeleton_items_2
{
public:

  template<class Refs, class Traits>
  struct Vertex_wrapper
  {
    typedef typename Traits::RT RT ;
    typedef typename Traits::Point_2 Point ;
    typedef Straight_skeleton_vertex_base_2 < Refs, Point, RT > Vertex; 
  };
  
  template<class Refs, class Traits> 
  struct Halfedge_wrapper
  {
    typedef typename Traits::Segment_2 Segment ;
    typedef Straight_skeleton_halfedge_base_2 < Refs, Segment > Halfedge; 
  };
  
  template<class Refs, class Traits> 
  struct Face_wrapper 
  {
    typedef HalfedgeDS_face_base< Refs > Face;
  } ;
};

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_ITEMS_2_H //
// EOF //
 
