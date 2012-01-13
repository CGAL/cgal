// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_BASE_FUNCTOR_H
#define CGAL_GPS_BASE_FUNCTOR_H

namespace CGAL {


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
  


   void create_face (Face_const_handle ,
                     Face_const_handle,
                     Face_handle )
  {}

  void create_vertex(Halfedge_const_handle ,
                     Halfedge_const_handle ,
                     Vertex_handle )
  {}

  void create_vertex(Vertex_const_handle ,
                     Vertex_const_handle ,
                     Vertex_handle  )
  {}

  void create_vertex(Vertex_const_handle ,
                     Halfedge_const_handle ,
                     Vertex_handle )
  {}

  void create_vertex(Halfedge_const_handle ,
                     Vertex_const_handle ,
                     Vertex_handle )
  {}

  void create_vertex(Face_const_handle ,
                     Vertex_const_handle ,
                     Vertex_handle )
  {}

  void create_vertex(Vertex_const_handle ,
                     Face_const_handle ,
                     Vertex_handle )
  {}

  void create_edge(Halfedge_const_handle ,
                   Halfedge_const_handle ,
                   Halfedge_handle )
  {}

  void create_edge(Halfedge_const_handle ,
                   Face_const_handle ,
                   Halfedge_handle )
  {}

  void create_edge(Face_const_handle ,
                   Halfedge_const_handle ,
                   Halfedge_handle )
  {}

};


} //namespace CGAL

#endif
