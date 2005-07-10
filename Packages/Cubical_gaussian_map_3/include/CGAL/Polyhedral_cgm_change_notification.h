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
// $Revision$
// $Name$
//
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_POLYHEDRAL_CGM_CHANGE_NOTOFICATION_H
#define CGAL_POLYHEDRAL_CGM_CHANGE_NOTOFICATION_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template<class T_Cgm>
class Polyhedral_cgm_change_notification {
public:
  typedef T_Cgm                                         Cgm;
  typedef typename Cgm::Polyhedron_vertex_const_handle
    Polyhedron_vertex_const_handle;
  typedef typename Cgm::Polyhedron_halfedge_const_handle
    Polyhedron_halfedge_const_handle;
  typedef typename Cgm::Polyhedron_facet_const_handle
    Polyhedron_facet_const_handle;

  typedef typename Cgm::Arr_vertex_handle               Arr_vertex_handle;
  typedef typename Cgm::Arr_halfedge_handle             Arr_halfedge_handle;
  typedef typename Cgm::Arr_face_handle                 Arr_face_handle;
  typedef typename Cgm::Arr_face_const_handle           Arr_face_const_handle;

  /*! Pass information from a polyhedron vertex to its dual - a cgm-face */
  virtual void update_dual_vertex(Polyhedron_vertex_const_handle src,
                                  Arr_face_handle trg)
  {}

  /*! Pass information from a dual vertex (face) to another dual vertex (face)
   * of the same vertex
   */
  virtual void update_dual_vertex(Arr_face_const_handle src,
                                  Arr_face_handle trg)
  {}

  /*! Pass information from a polyhedron facet to its dual a cgm-vertex */
  virtual void update_dual_face(Polyhedron_facet_const_handle src,
                                Arr_vertex_handle trg)
  {}

  /*! Pass information from a polyhedron facet to its dual a cgm-vertex */
  virtual void update_dual_halfedge(Polyhedron_halfedge_const_handle src,
                                    Arr_halfedge_handle trg)
  {}
};

CGAL_END_NAMESPACE

#endif
