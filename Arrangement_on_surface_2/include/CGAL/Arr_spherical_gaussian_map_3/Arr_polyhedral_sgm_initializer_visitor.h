// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POLYHEDRAL_SGM_INITIALIZER_VISITOR_H
#define CGAL_ARR_POLYHEDRAL_SGM_INITIALIZER_VISITOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/basic.h>

namespace CGAL {

template<class PolyhedralSgm, class Polyhedron>
class Arr_polyhedral_sgm_initializer_visitor {
public:
  typedef typename Polyhedron::Vertex_const_handle
                                        Polyhedron_vertex_const_handle;
  typedef typename Polyhedron::Halfedge_const_handle
                                        Polyhedron_halfedge_const_handle;
  typedef typename Polyhedron::Facet_const_handle
                                        Polyhedron_facet_const_handle;

  typedef typename PolyhedralSgm::Vertex_handle         Sgm_vertex_handle;
  typedef typename PolyhedralSgm::Halfedge_handle       Sgm_halfedge_handle;
  typedef typename PolyhedralSgm::Face_handle           Sgm_face_handle;
  typedef typename PolyhedralSgm::Face_const_handle     Sgm_face_const_handle;

  /*! Destructor */
  virtual ~Arr_polyhedral_sgm_initializer_visitor() {}

  /*! Pass information from a polyhedron vertex to its dual - a sgm-face */
  virtual void update_dual_vertex(Polyhedron_vertex_const_handle /*src*/,
                                  Sgm_face_handle /*trg*/)
  {}

  /*! Pass information from a dual vertex (face) to another dual vertex (face)
   * of the same vertex
   */
  virtual void update_dual_vertex(Sgm_face_const_handle /*src*/,
                                  Sgm_face_handle /*trg*/)
  {}

  /*! Pass information from a polyhedron facet to its dual a sgm-vertex */
  virtual void update_dual_face(Polyhedron_facet_const_handle /*src*/,
                                Sgm_vertex_handle /*trg*/)
  {}

  /*! Pass information from a polyhedron facet to its dual a sgm-vertex */
  virtual void update_dual_halfedge(Polyhedron_halfedge_const_handle /*src*/,
                                    Sgm_halfedge_handle /*trg*/)
  {}
};

} //namespace CGAL

#endif
