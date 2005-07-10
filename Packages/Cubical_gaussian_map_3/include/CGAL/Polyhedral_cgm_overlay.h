// Copyright(c) 2005  Tel-Aviv University(Israel).
// All rights reserved.
//
// This file is part of CGAL(www.cgal.org); you may redistribute it under
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

#ifndef CGAL_POLYHEDRAL_CGM_OVERLAY_H
#define CGAL_POLYHEDRAL_CGM_OVERLAY_H

#include <iostream>

CGAL_BEGIN_NAMESPACE

template <class T_Cgm>
class Polyhedral_cgm_overlay {
private:
  typedef T_Cgm                                         Cgm;
  typedef typename Cgm::Kernel                          Kernel;
  typedef typename Cgm::Arrangement                     Arrangement;
  typedef typename Cgm::Point_3                         Point_3;
  typedef typename Cgm::Vector_3                        Vector_3;
  typedef typename Cgm::Arr_vertex                      Arr_vertex;
  
  unsigned int m_face_id;
  
public:
  /*! Obtain the index of the containing unit-cube face */
  unsigned int get_face_id() const { return m_face_id; }

  /*! Set the index of the containing unit-cube face */
  void set_face_id(unsigned int face_id) { m_face_id = face_id; }
  
  typedef typename Arrangement::Face_handle             Face_handle;
  typedef typename Arrangement::Vertex_handle           Vertex_handle;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;

  typedef typename Arrangement::Face_const_handle       Face_const_handle;
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;

  typedef typename Arrangement::Ccb_halfedge_const_circulator 
    Ccb_halfedge_const_circulator;

  typedef typename Arrangement::Halfedge_around_vertex_const_circulator
    Arr_halfedge_around_vertex_const_circulator;
  
  /*! 1 */
  void create_face(Face_const_handle f1, Face_const_handle f2,
                   Face_handle f)
  {
    const Point_3 & p1 = f1->get_point();
    const Point_3 & p2 = f2->get_point();
    Vector_3 v1(ORIGIN, p1);
    Point_3 p = p2 + v1;
    f->set_point(p);
  }

  /*! 2 */
  void create_vertex(Halfedge_const_handle h1, Halfedge_const_handle h2,
                     Vertex_handle v)
  {
    v->set_is_real(true);
    v->set_face_id(m_face_id);
    v->set_location(Arr_vertex::Interior);
  }

  /*! 3 */
  void create_vertex(Vertex_const_handle v1, Vertex_const_handle v2,
                     Vertex_handle v)
  {
    v->set_face_id(m_face_id);
    v->set_location(v1->get_location());

    // If one of the vertices is real, the overlay vertex is real, or
    v->set_is_real(v1->get_is_real() || v2->get_is_real());

    if (v1->get_location() == Arr_vertex::Interior) return;

    /* If the vertex is a Corner or an Edge vertex or a Corner vertex, then it
     * is real iff both input vertices have real incident halfedges that do not
     * overlap
     */
    bool is_real1 = false;
    Arr_halfedge_around_vertex_const_circulator hec1 =
      v1->incident_halfedges();
    Arr_halfedge_around_vertex_const_circulator begin_hec1 = hec1;
    do {
      if (hec1->get_is_real()) {
        is_real1 = true;
        break;
      }
      ++hec1;
    } while (hec1 != begin_hec1);

    bool is_real2 = false;
    Arr_halfedge_around_vertex_const_circulator hec2 =
      v2->incident_halfedges();
    Arr_halfedge_around_vertex_const_circulator begin_hec2 = hec2;
    do {
      if (hec2->get_is_real()) {
        is_real2 = true;
        break;
      }
      ++hec2;
    } while (hec2 != begin_hec2);

    if (is_real1 && is_real2) {
      Kernel kernel;
      if (kernel.orientation_2_object()(hec1->source()->point(),
                                        hec2->source()->point(),
                                        v->point()) != COLLINEAR)
        v->set_is_real(true);
    }
  }

  /*! 4 */
  void create_vertex(Vertex_const_handle v1, Halfedge_const_handle h2,
                     Vertex_handle v)
  {
    v->set_face_id(m_face_id);
    v->set_is_real(v1->get_is_real() || h2->get_is_real());
    v->set_location(v1->get_location());
  }

  /*! 5 */
  void create_vertex(Halfedge_const_handle h1, Vertex_const_handle v2,
                     Vertex_handle v)
  {
    v->set_face_id(m_face_id);
    v->set_is_real(v2->get_is_real() || h1->get_is_real());
    v->set_location(v2->get_location());
  }

  /*! 6 */
  void create_vertex(Face_const_handle f1, Vertex_const_handle v2,
                     Vertex_handle v)
  {
    v->set_face_id(m_face_id);
    v->set_is_real(v2->get_is_real());
    v->set_location(v2->get_location());
  }

  /*! 7 */
  void create_vertex(Vertex_const_handle v1, Face_const_handle f2,
                     Vertex_handle v)
  {
    v->set_face_id(m_face_id);
    v->set_is_real(v1->get_is_real());
    v->set_location(v1->get_location());
  }

  /*! 8 */
  void create_edge(Halfedge_const_handle h1, Halfedge_const_handle h2,
                   Halfedge_handle h)
  {
    bool is_real = h1->get_is_real() || h2->get_is_real();
    h->set_is_real(is_real);
    h->twin()->set_is_real(is_real);
  }

  /*! 9 */
  void create_edge(Halfedge_const_handle h1, Face_const_handle f2,
                   Halfedge_handle h)
  {
    bool is_real = h1->get_is_real();
    h->set_is_real(is_real);
    h->twin()->set_is_real(is_real);
  }

  /*! 10 */
  void create_edge(Face_const_handle f1, Halfedge_const_handle h2,
                   Halfedge_handle h)
  {
    bool is_real = h2->get_is_real();
    h->set_is_real(is_real);
    h->twin()->set_is_real(is_real);
  }
};

CGAL_END_NAMESPACE

#endif
