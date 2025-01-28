// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_DO_INTERSECT_FUNCTOR_H
#define CGAL_GPS_DO_INTERSECT_FUNCTOR_H

#include <CGAL/license/Boolean_set_operations_2.h>


namespace CGAL {

template <class Arrangement_>
class Gps_do_intersect_functor
{
public:

  typedef Arrangement_                                    Arrangement_2;

  typedef typename Arrangement_2::Face_const_handle       Face_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle   Halfedge_const_handle;

  typedef typename Arrangement_2::Face_handle        Face_handle;
  typedef typename Arrangement_2::Halfedge_handle    Halfedge_handle;
  typedef typename Arrangement_2::Vertex_handle      Vertex_handle;

  // default constructor
  Gps_do_intersect_functor() : m_found_reg_intersection(false),
                               m_found_boudary_intersection(false)

  {}

   void create_face (Face_const_handle f1,
                     Face_const_handle f2,
                     Face_handle )
  {
    if(f1->contained() && f2->contained())
      // found intersection
      m_found_reg_intersection = true;
  }


  void create_vertex(Vertex_const_handle ,
                     Vertex_const_handle ,
                     Vertex_handle  )
  {
    m_found_boudary_intersection = true;
  }

  void create_vertex(Vertex_const_handle ,
                     Halfedge_const_handle ,
                     Vertex_handle )
  {
    m_found_boudary_intersection = true;
  }

  void create_vertex(Halfedge_const_handle ,
                     Vertex_const_handle ,
                     Vertex_handle )
  {
    m_found_boudary_intersection = true;
  }

  void create_vertex(Halfedge_const_handle ,
                     Halfedge_const_handle ,
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
  {
    m_found_boudary_intersection = true;
  }

  void create_edge(Halfedge_const_handle ,
                   Face_const_handle ,
                   Halfedge_handle )
  {}

  void create_edge(Face_const_handle ,
                   Halfedge_const_handle ,
                   Halfedge_handle )
  {}


  bool found_reg_intersection() const
  {
    return m_found_reg_intersection;
  }

  bool found_boundary_intersection() const
  {
    return m_found_boudary_intersection;
  }

  protected:

    bool m_found_reg_intersection;
    bool m_found_boudary_intersection;
};


} //namespace CGAL

#endif
