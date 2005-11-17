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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef BSO_DO_INTERSECT_FUNCTOR
#define BSO_DO_INTERSECT_FUNCTOR

#include <CGAL/Arrangement_2.h>
#include <CGAL/Boolean_set_operations_2/Bso_dcel.h>

CGAL_BEGIN_NAMESPACE

template <class Traits_>
class Bso_do_intersect_functor
{
public:

  typedef Traits_                        Traits;
  typedef Bso_dcel<Traits>              Dcel;
  typedef Arrangement_2<Traits, Dcel>    Bso_arrangement;

  typedef typename Bso_arrangement::Face_const_handle       Face_const_handle;
  typedef typename Bso_arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Bso_arrangement::Halfedge_const_handle   Halfedge_const_handle;

  typedef typename Bso_arrangement::Face_handle        Face_handle;
  typedef typename Bso_arrangement::Halfedge_handle    Halfedge_handle;
  typedef typename Bso_arrangement::Vertex_handle      Vertex_handle;


  // default constructor
  Bso_do_intersect_functor() : m_found_x(false)
  {}

   void create_face (Face_const_handle f1,
                     Face_const_handle f2,
                     Face_handle res_f)
  {
    if(f1->contained() && f2->contained())
      // found intersection
      m_found_x = true;    
  }

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

  bool found_intersection() const
  {
    return m_found_x;
  }

  private:

    bool m_found_x;
};


CGAL_END_NAMESPACE

#endif
