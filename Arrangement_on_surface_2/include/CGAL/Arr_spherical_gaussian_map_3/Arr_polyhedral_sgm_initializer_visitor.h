// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
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
