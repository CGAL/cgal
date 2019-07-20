// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_FACE_BASE_2_H
#define CGAL_STRAIGHT_SKELETON_FACE_BASE_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/tags.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL {

template <class Refs>
class Straight_skeleton_face_base_base_2
{
public:

  typedef Straight_skeleton_face_base_base_2<Refs> Base;
  
  typedef Refs                                 HalfedgeDS;
  typedef Tag_true                             Supports_face_halfedge;
  typedef typename Refs::Vertex_handle         Vertex_handle;
  typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
  typedef typename Refs::Halfedge_handle       Halfedge_handle;
  typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Refs::Face_handle           Face_handle;
  typedef typename Refs::Face_const_handle     Face_const_handle;
  typedef typename Refs::Vertex                Vertex;
  typedef typename Refs::Halfedge              Halfedge;

  Straight_skeleton_face_base_base_2() : mID(-1) {}

  Straight_skeleton_face_base_base_2( int aID ) : mID(aID) {}
  
public:

  int id() const { return mID ; }
  
  Halfedge_handle       halfedge()       { return mHE; }
  Halfedge_const_handle halfedge() const { return mHE; }
  
  void set_halfedge( Halfedge_handle aHE )  { mHE = aHE; }
  
  void reset_id ( int aID ) { mID = aID ; }
  
private:
  
  int             mID ;
  Halfedge_handle mHE;
  
};

template < class Refs >
class Straight_skeleton_face_base_2 : public Straight_skeleton_face_base_base_2<Refs>
{
public:

  typedef typename Refs::Vertex_handle   Vertex_handle;
  typedef typename Refs::Halfedge_handle Halfedge_handle;
  typedef typename Refs::Face_handle     Face_handle;
  
  typedef Straight_skeleton_face_base_base_2<Refs> Base ;
  
  Straight_skeleton_face_base_2() {}
  
  Straight_skeleton_face_base_2( int aID ) : Base(aID) {}

private:

  void set_halfedge( Halfedge_handle aHE ) { Base::set_halfedge(aHE) ; }
  void reset_id    ( int aID )             { Base::reset_id(aID) ; }

} ;

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_FACE_BASE_2_H //
// EOF //

