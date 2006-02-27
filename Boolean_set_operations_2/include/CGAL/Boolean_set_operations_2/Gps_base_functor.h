// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_BASE_FUNCTOR_H
#define CGAL_GPS_BASE_FUNCTOR_H

CGAL_BEGIN_NAMESPACE


template <class Arrangement_>
class Gps_base_functor
{
public:

  typedef Arrangement_       Arrangement_2;

  typedef typename Arrangement_2::Face_const_handle       Face_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle   Halfedge_const_handle;

  typedef typename Arrangement_2::Face_handle        Face_handle;
  typedef typename Arrangement_2::Halfedge_handle    Halfedge_handle;
  typedef typename Arrangement_2::Vertex_handle      Vertex_handle;
  


   void create_face (Face_const_handle f1,
                     Face_const_handle f2,
                     Face_handle res_f)
  {}

  void create_vertex(Halfedge_const_handle h1,
                     Halfedge_const_handle h2,
                     Vertex_handle res_v)
  {}

  void create_vertex(Vertex_const_handle v1,
                     Vertex_const_handle v2,
                     Vertex_handle  res_v)
  {}

  void create_vertex(Vertex_const_handle v1,
                     Halfedge_const_handle h2,
                     Vertex_handle res_v)
  {}

  void create_vertex(Halfedge_const_handle h1,
                     Vertex_const_handle v2,
                     Vertex_handle res_v)
  {}

  void create_vertex(Face_const_handle f1,
                     Vertex_const_handle v2,
                     Vertex_handle res_v)
  {}

  void create_vertex(Vertex_const_handle v1,
                     Face_const_handle f2,
                     Vertex_handle res_v)
  {}

  void create_edge(Halfedge_const_handle h1,
                   Halfedge_const_handle h2,
                   Halfedge_handle res_h)
  {}

  void create_edge(Halfedge_const_handle h1,
                   Face_const_handle f2,
                   Halfedge_handle res_h)
  {}

  void create_edge(Face_const_handle f1,
                   Halfedge_const_handle h2,
                   Halfedge_handle res_h)
  {}

};


CGAL_END_NAMESPACE

#endif
