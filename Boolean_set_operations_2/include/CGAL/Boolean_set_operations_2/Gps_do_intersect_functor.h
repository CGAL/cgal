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

template <typename Arrangement_>
class Gps_do_intersect_functor {
public:
  using Arrangement_2 = Arrangement_;

  using Face_const_handle = typename Arrangement_2::Face_const_handle;
  using Vertex_const_handle = typename Arrangement_2::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement_2::Halfedge_const_handle;

  using Face_handle = typename Arrangement_2::Face_handle;
  using Halfedge_handle = typename Arrangement_2::Halfedge_handle;
  using Vertex_handle = typename Arrangement_2::Vertex_handle;

  // default constructor
  Gps_do_intersect_functor() :
    m_found_reg_intersection(false),
    m_found_boudary_intersection(false)
  {}

  void create_face(Face_const_handle f1, Face_const_handle f2, Face_handle)
  { if (f1->contained() && f2->contained()) m_found_reg_intersection = true; }

  void create_vertex(Vertex_const_handle, Vertex_const_handle, Vertex_handle)
  { m_found_boudary_intersection = true; }

  void create_vertex(Vertex_const_handle, Halfedge_const_handle, Vertex_handle)
  { m_found_boudary_intersection = true; }

  void create_vertex(Halfedge_const_handle, Vertex_const_handle, Vertex_handle)
  { m_found_boudary_intersection = true; }

  void create_vertex(Halfedge_const_handle, Halfedge_const_handle, Vertex_handle) {}

  void create_vertex(Face_const_handle, Vertex_const_handle, Vertex_handle) {}

  void create_vertex(Vertex_const_handle, Face_const_handle, Vertex_handle) {}

  void create_edge(Halfedge_const_handle, Halfedge_const_handle, Halfedge_handle)
  { m_found_boudary_intersection = true; }

  void create_edge(Halfedge_const_handle, Face_const_handle, Halfedge_handle) {}

  void create_edge(Face_const_handle, Halfedge_const_handle, Halfedge_handle) {}

  bool found_reg_intersection() const { return m_found_reg_intersection; }

  bool found_boundary_intersection() const { return m_found_boudary_intersection; }

protected:
  bool m_found_reg_intersection;
  bool m_found_boudary_intersection;
};

} //namespace CGAL

#endif
