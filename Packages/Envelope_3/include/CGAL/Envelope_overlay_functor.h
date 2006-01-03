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
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef ENVELOPE_OVERLAY_FUNCTOR_H
#define ENVELOPE_OVERLAY_FUNCTOR_H

#include <iostream>

CGAL_BEGIN_NAMESPACE

template < class MinimizationDiagram_2>
class Envelope_overlay_functor
{
public:
  typedef MinimizationDiagram_2                                  Minimization_diagram_2;
  
  typedef typename Minimization_diagram_2::Face_const_handle     Face_const_handle1;
  typedef typename Minimization_diagram_2::Face_const_handle     Face_const_handle2;

  typedef typename Minimization_diagram_2::Vertex_const_handle   Vertex_const_handle1;
  typedef typename Minimization_diagram_2::Vertex_const_handle   Vertex_const_handle2;

  typedef typename Minimization_diagram_2::Halfedge_const_handle Halfedge_const_handle1;
  typedef typename Minimization_diagram_2::Halfedge_const_handle Halfedge_const_handle2;

  typedef typename Minimization_diagram_2::Face_handle           Face_handle;
  typedef typename Minimization_diagram_2::Vertex_handle         Vertex_handle;
  typedef typename Minimization_diagram_2::Halfedge_handle       Halfedge_handle;

  typedef typename Minimization_diagram_2::Face_handle           Res_face_handle;
  typedef typename Minimization_diagram_2::Halfedge_handle       Res_halfedge_handle;
  typedef typename Minimization_diagram_2::Vertex_handle         Res_vertex_handle;

protected:

  typedef typename Minimization_diagram_2::Dcel::Dcel_data_iterator Envelope_data_iterator;
public:

  Envelope_overlay_functor(Minimization_diagram_2& md1,
                           Minimization_diagram_2& md2,
                           Minimization_diagram_2& result)
    : m_1(md1), m_2(md2), m_result(result)
  {
  }
  
  void create_face (Face_const_handle1 f1, Face_const_handle2 f2, Res_face_handle res_f)
  {
    res_f->set_aux_source(0, m_1.non_const_handle(f1));
    res_f->set_aux_source(1, m_2.non_const_handle(f2));
  }

  void create_vertex(Halfedge_const_handle1 h1,
                     Halfedge_const_handle2 h2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(h1));
    res_v->set_aux_source(1, m_2.non_const_handle(h2));
    res_v->set_is_intersection(true);
    // res_v cannot be isolated
  }

  void create_vertex(Vertex_const_handle1 v1,
                     Vertex_const_handle2 v2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(v1));
    res_v->set_aux_source(1, m_2.non_const_handle(v2));
    res_v->set_is_intersection(false);

    if (v1->is_isolated() && v2->is_isolated())
    {
      res_v->set_is_equal_aux_data_in_face(0, v1->get_is_equal_data_in_face());
      res_v->set_is_equal_aux_data_in_face(1, v2->get_is_equal_data_in_face());
      res_v->set_has_equal_aux_data_in_face(0, v1->get_has_equal_data_in_face());
      res_v->set_has_equal_aux_data_in_face(1, v2->get_has_equal_data_in_face());      
    }
  }

  void create_vertex(Vertex_const_handle1 v1,
                     Halfedge_const_handle2 h2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(v1));
    res_v->set_aux_source(1, m_2.non_const_handle(h2));
    res_v->set_is_intersection(true);
    // res_v cannot be isolated
  }

  void create_vertex(Halfedge_const_handle1 h1,
                     Vertex_const_handle2 v2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(h1));
    res_v->set_aux_source(1, m_2.non_const_handle(v2));
    res_v->set_is_intersection(true);
    // res_v cannot be isolated
  }

  void create_vertex(Face_const_handle1 f1,
                     Vertex_const_handle2 v2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(f1));
    res_v->set_aux_source(1, m_2.non_const_handle(v2));
    res_v->set_is_intersection(false);

    if (v2->is_isolated())
    {
      // the res_v is also isolated, and we should update the is_equal/has_equal
      // data in face information
      res_v->set_is_equal_aux_data_in_face(0, true);
      res_v->set_is_equal_aux_data_in_face(1, v2->get_is_equal_data_in_face());
      res_v->set_has_equal_aux_data_in_face(0, !f1->has_no_data());
      res_v->set_has_equal_aux_data_in_face(1, v2->get_has_equal_data_in_face());
    }
  }

  void create_vertex(Vertex_const_handle1 v1,
                     Face_const_handle2 f2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(v1));
    res_v->set_aux_source(1, m_2.non_const_handle(f2));
    res_v->set_is_intersection(false);

    if (v1->is_isolated())
    {
      // the res_v is also isolated, and we should update the is_equal/has_equal
      // data in face information
      res_v->set_is_equal_aux_data_in_face(0, v1->get_is_equal_data_in_face());
      res_v->set_is_equal_aux_data_in_face(1, true);
      res_v->set_has_equal_aux_data_in_face(0, v1->get_has_equal_data_in_face());
      res_v->set_has_equal_aux_data_in_face(1, !f2->has_no_data());      
    }
  }

  void create_edge(Halfedge_const_handle1 h1,
                   Halfedge_const_handle2 h2,
                   Res_halfedge_handle res_h)
  {
    // update source
    res_h->set_aux_source(0, m_1.non_const_handle(h1));
    res_h->set_aux_source(1, m_2.non_const_handle(h2));

    res_h->twin()->set_aux_source(0, m_1.non_const_handle(h1->twin()));
    res_h->twin()->set_aux_source(1, m_2.non_const_handle(h2->twin()));

    // update is_equal/has_equal data in face
    res_h->set_is_equal_aux_data_in_face(0, h1->get_is_equal_data_in_face());
    res_h->set_is_equal_aux_data_in_face(1, h2->get_is_equal_data_in_face());
    res_h->set_has_equal_aux_data_in_face(0, h1->get_has_equal_data_in_face());
    res_h->set_has_equal_aux_data_in_face(1, h2->get_has_equal_data_in_face());

    res_h->twin()->set_is_equal_aux_data_in_face(0, h1->twin()->get_is_equal_data_in_face());
    res_h->twin()->set_is_equal_aux_data_in_face(1, h2->twin()->get_is_equal_data_in_face());
    res_h->twin()->set_has_equal_aux_data_in_face(0, h1->twin()->get_has_equal_data_in_face());
    res_h->twin()->set_has_equal_aux_data_in_face(1, h2->twin()->get_has_equal_data_in_face());

    // update is_equal/has_equal data in target
    Vertex_handle vh1;
    Halfedge_handle hh1;
    // update aux_data(0)
    Object trg_src1 = res_h->target()->get_aux_source(0);
    if (assign(vh1, trg_src1))
      // v is the target of h1, and we can copy the halfedge-target information
      // from h1
      copy_halfedge_target_info(h1, res_h, 0);
    else if (assign(hh1, trg_src1))
      // h is the "HEMSHECH" of h1, so we need to set halfedge_target
      // information to true
      set_halfedge_target_info(res_h, 0, true);
    else
      // this cannot happen, since we need to touch an edge
      CGAL_assertion(false);

    // update aux_data(1)
    Vertex_handle vh2;
    Halfedge_handle hh2;
    Object trg_src2 = res_h->target()->get_aux_source(1);
    if (assign(vh2, trg_src2))
      // v is the target of h2, and we can copy the halfedge-target information
      // from h1
      copy_halfedge_target_info(h2, res_h, 1);
    else if (assign(hh2, trg_src2))
      // h is the "HEMSHECH" of h2, so we need to set halfedge_target
      // information to true
      set_halfedge_target_info(res_h, 1, true);
    else
      // this cannot happen, since we need to touch an edge
      CGAL_assertion(false);

    // update is_equal/has_equal data in source
    // update aux_data(0)
    Object src_src1 = res_h->source()->get_aux_source(0);
    if (assign(vh1, src_src1))
      // v is the target of h1->twin(), and we can copy the halfedge-target
      // information from h1->twin()
      copy_halfedge_target_info(h1->twin(), res_h->twin(), 0);
    else if (assign(hh1, src_src1))
      // h is the "HEMSHECH" of h1->twin(), so we need to set halfedge_target
      // information to true
      set_halfedge_target_info(res_h->twin(), 0, true);
    else
      // this cannot happen, since we need to touch an edge
      CGAL_assertion(false);

    // update aux_data(1)
    Object src_src2 = res_h->source()->get_aux_source(1);
    if (assign(vh2, src_src2))
      // v is the target of h2->twin(), and we can copy the halfedge-target
      // information from h2->twin()
      copy_halfedge_target_info(h2->twin(), res_h->twin(), 1);
    else if (assign(hh2, src_src2))
      // h is the "HEMSHECH" of h2->twin(), so we need to set halfedge_target
      // information to true
      set_halfedge_target_info(res_h->twin(), 1, true);
    else
      // this cannot happen, since we need to touch an edge
      CGAL_assertion(false);

  }

  void create_edge(Halfedge_const_handle1 h1,
                   Face_const_handle2 f2,
                   Res_halfedge_handle res_h)
  {
    res_h->set_aux_source(0, m_1.non_const_handle(h1));
    res_h->set_aux_source(1, m_2.non_const_handle(f2));

    res_h->twin()->set_aux_source(0, m_1.non_const_handle(h1->twin()));
    res_h->twin()->set_aux_source(1, m_2.non_const_handle(f2));

    res_h->set_is_equal_aux_data_in_face(0, h1->get_is_equal_data_in_face());
    res_h->set_is_equal_aux_data_in_face(1, true);
    res_h->set_has_equal_aux_data_in_face(0, h1->get_has_equal_data_in_face());
    res_h->set_has_equal_aux_data_in_face(1, !f2->has_no_data());

    res_h->twin()->set_is_equal_aux_data_in_face(0, h1->twin()->get_is_equal_data_in_face());
    res_h->twin()->set_is_equal_aux_data_in_face(1, true);
    res_h->twin()->set_has_equal_aux_data_in_face(0, h1->twin()->get_has_equal_data_in_face());
    res_h->twin()->set_has_equal_aux_data_in_face(1, !f2->has_no_data());

    // update is_equal/has_equal data in target for the fisrt source map 
    Vertex_handle v;
    Halfedge_handle h;
    // update target
    Object trg_src1 = res_h->target()->get_aux_source(0);
    if (assign(v, trg_src1))
      // v is the target of h1, and we can copy the halfedge-target information
      // from h1
      copy_halfedge_target_info(h1, res_h, 0);
    else if (assign(h, trg_src1))
      // h is the "HEMSHECH" of h1, so we need to set halfedge_target
      // information to true
      set_halfedge_target_info(res_h, 0, true);
    else
      // this cannot happen, since we need to touch an edge
      CGAL_assertion(false);

    // update source
    Object src_src1 = res_h->source()->get_aux_source(0);
    if (assign(v, src_src1))
      // v is the target of h1->twin(), and we can copy the halfedge-target
      // information from h1->twin()
      copy_halfedge_target_info(h1->twin(), res_h->twin(), 0);
    else if (assign(h, src_src1))
      // h is the "HEMSHECH" of h1->twin(), so we need to set halfedge_target
      // information to true
      set_halfedge_target_info(res_h->twin(), 0, true);
    else
      // this cannot happen, since we need to touch an edge
      CGAL_assertion(false);

    // update is_equal/has_equal data in target for the second source map
    Vertex_handle vh2;
    Halfedge_handle hh2;
    Face_handle fh2;
    // update target
    Object trg_src2 = res_h->target()->get_aux_source(1);
    if (assign(vh2, trg_src2))
      // todo
    {
      if (vh2->is_isolated())
        copy_halfedge_target_info_from_vertex_face_info(vh2, res_h, 1);
      else
      {
        // we have a vertex vh2 on the boundary of the face f2
        bool is_equal = vh2->is_equal_data(f2->begin_data(), f2->end_data());
        bool has_equal = vh2->has_equal_data(f2->begin_data(), f2->end_data());

        res_h->set_is_equal_aux_data_in_target(1, is_equal);
        res_h->set_has_equal_aux_data_in_target(1, has_equal);
      }
    } 
    else if (assign(hh2, trg_src2))
      // we should find the halfedge (hh2 or hh2->twin()) that points to face
      // f2, and check the in_face flags there
    {
      CGAL_assertion(hh2->face() == m_2.non_const_handle(f2) ||
                     hh2->twin()->face() == m_2.non_const_handle(f2));
      Face_handle f = m_2.non_const_handle(f2);
      if (hh2->face() == f)
        copy_halfedge_target_info_from_halfedge_face_info(hh2, res_h, 1);
      else
        copy_halfedge_target_info_from_halfedge_face_info(hh2->twin(), res_h, 1);
    }
    else
    {
      CGAL_assertion(assign(fh2, trg_src2));
      assign(fh2, trg_src2);
      // the edge and target are in the same face, so we set halfedge-target
      // equal information to true, and has equal acoording to face data
      res_h->set_is_equal_aux_data_in_target(1, true);
      res_h->set_has_equal_aux_data_in_target(1, !fh2->has_no_data());
//      set_halfedge_target_info(res_h, 1, true);
    }
    
    // update source
    Object src_src2 = res_h->source()->get_aux_source(1);
    if (assign(vh2, src_src2))
      // todo
    {
      // we have a vertex vh2 on the boundary of the face f2
      // or an isolated vertex
      if (vh2->is_isolated())
        copy_halfedge_target_info_from_vertex_face_info(vh2, res_h->twin(), 1);
      else
      {
        bool is_equal = vh2->is_equal_data(f2->begin_data(), f2->end_data());
        bool has_equal = vh2->has_equal_data(f2->begin_data(), f2->end_data());

        res_h->twin()->set_is_equal_aux_data_in_target(1, is_equal);
        res_h->twin()->set_has_equal_aux_data_in_target(1, has_equal);
      }
    }
    else if (assign(hh2, src_src2))
      // we should find the halfedge (hh2 or hh2->twin()) that points to face
      // f2, and check the in_face flags there
    {
      CGAL_assertion(hh2->face() == m_2.non_const_handle(f2) ||
                     hh2->twin()->face() == m_2.non_const_handle(f2));
      Face_handle f = m_2.non_const_handle(f2);
      if (hh2->face() == f)
        copy_halfedge_target_info_from_halfedge_face_info(hh2, res_h->twin(), 1);
      else
        copy_halfedge_target_info_from_halfedge_face_info(hh2->twin(), res_h->twin(), 1);
    }

    else
    {
      CGAL_assertion(assign(fh2, src_src2));
      assign(fh2, src_src2);
      // the edge and soucre are in the same face, so we set halfedge-target
      // information to true
      res_h->twin()->set_is_equal_aux_data_in_target(1, true);
      res_h->twin()->set_has_equal_aux_data_in_target(1, !fh2->has_no_data());

//      set_halfedge_target_info(res_h->twin(), 1, true);
    }

  }                           


  void create_edge(Face_const_handle1 f1,
                   Halfedge_const_handle2 h2,
                   Res_halfedge_handle res_h)
  {
    res_h->set_aux_source(0, m_1.non_const_handle(f1));
    res_h->set_aux_source(1, m_2.non_const_handle(h2));

    res_h->twin()->set_aux_source(0, m_1.non_const_handle(f1));
    res_h->twin()->set_aux_source(1, m_2.non_const_handle(h2->twin()));

    res_h->set_is_equal_aux_data_in_face(0, true);
    res_h->set_is_equal_aux_data_in_face(1, h2->get_is_equal_data_in_face());
    res_h->set_has_equal_aux_data_in_face(0, !f1->has_no_data());
    res_h->set_has_equal_aux_data_in_face(1, h2->get_has_equal_data_in_face());

    res_h->twin()->set_is_equal_aux_data_in_face(0, true);
    res_h->twin()->set_is_equal_aux_data_in_face(1, h2->twin()->get_is_equal_data_in_face());
    res_h->twin()->set_has_equal_aux_data_in_face(0, !f1->has_no_data());
    res_h->twin()->set_has_equal_aux_data_in_face(1, h2->twin()->get_has_equal_data_in_face());

    // update is_equal/has_equal data in target for the second source map
    Vertex_handle v;
    Halfedge_handle h;
    // update target
    Object trg_src2 = res_h->target()->get_aux_source(1);
    if (assign(v, trg_src2))
      // v is the target of h2, and we can copy the halfedge-target information
      // from h1
      copy_halfedge_target_info(h2, res_h, 1);
    else if (assign(h, trg_src2))
      // h is the "HEMSHECH" of h2, so we need to set halfedge_target
      // information to true
      set_halfedge_target_info(res_h, 1, true);
    else
      // this cannot happen, since we need to touch an edge
      CGAL_assertion(false);

    // update source
    Object src_src2 = res_h->source()->get_aux_source(1);
    if (assign(v, src_src2))
      // v is the target of h2->twin(), and we can copy the halfedge-target
      // information from h2->twin()
      copy_halfedge_target_info(h2->twin(), res_h->twin(), 1);
    else if (assign(h, src_src2))
      // h is the "HEMSHECH" of h2->twin(), so we need to set halfedge_target
      // information to true
      set_halfedge_target_info(res_h->twin(), 1, true);
    else
      // this cannot happen, since we need to touch an edge
      CGAL_assertion(false);

    // update is_equal/has_equal data in target for the first source map
    Vertex_handle vh1;
    Halfedge_handle hh1;
    Face_handle fh1;
    // update target
    Object trg_src1 = res_h->target()->get_aux_source(0);
    if (assign(vh1, trg_src1))
      // todo
    {
      if (vh1->is_isolated())
        copy_halfedge_target_info_from_vertex_face_info(vh1, res_h, 0);
      else
      {
        // we have a vertex vh1 on the boundary of the face f1
        bool is_equal = vh1->is_equal_data(f1->begin_data(), f1->end_data());
        bool has_equal = vh1->has_equal_data(f1->begin_data(), f1->end_data());

        res_h->set_is_equal_aux_data_in_target(0, is_equal);
        res_h->set_has_equal_aux_data_in_target(0, has_equal);
      }
    }
    else if (assign(hh1, trg_src1))
      // we should find the halfedge (hh1 or hh1>twin()) that points to face
      // f1, and check the in_face flags there
    {
      CGAL_assertion(hh1->face() == m_1.non_const_handle(f1) ||
                     hh1->twin()->face() == m_1.non_const_handle(f1));
      Face_handle f = m_1.non_const_handle(f1);
      if (hh1->face() == f)
        copy_halfedge_target_info_from_halfedge_face_info(hh1, res_h, 0);
      else
        copy_halfedge_target_info_from_halfedge_face_info(hh1->twin(), res_h, 0);
    }
    else
    {
      CGAL_assertion(assign(fh1, trg_src1));
      assign(fh1, trg_src1);
      // the edge and target are in the same face, so we set halfedge-target
      // information to true
      res_h->set_is_equal_aux_data_in_target(0, true);
      res_h->set_has_equal_aux_data_in_target(0, !fh1->has_no_data());
//      set_halfedge_target_info(res_h, 0, true);
    }

    // update source
    Object src_src1 = res_h->source()->get_aux_source(0);
    if (assign(vh1, src_src1))
      // todo
    {
      if (vh1->is_isolated())
        copy_halfedge_target_info_from_vertex_face_info(vh1, res_h->twin(), 0);
      else
      {
        // we have a vertex vh1 on the boundary of the face f1
        bool is_equal = vh1->is_equal_data(f1->begin_data(), f1->end_data());
        bool has_equal = vh1->has_equal_data(f1->begin_data(), f1->end_data());

        res_h->twin()->set_is_equal_aux_data_in_target(0, is_equal);
        res_h->twin()->set_has_equal_aux_data_in_target(0, has_equal);
      }
    }
    else if (assign(hh1, src_src1))
      // we should find the halfedge (hh1 or hh1->twin()) that points to face
      // f1, and check the in_face flags there
    {
      CGAL_assertion(hh1->face() == m_1.non_const_handle(f1) ||
                     hh1->twin()->face() == m_1.non_const_handle(f1));
      Face_handle f = m_1.non_const_handle(f1);
      if (hh1->face() == f)
        copy_halfedge_target_info_from_halfedge_face_info(hh1, res_h->twin(), 0);
      else
        copy_halfedge_target_info_from_halfedge_face_info(hh1->twin(), res_h->twin(), 0);
    }

    else
    {
      CGAL_assertion(assign(fh1, src_src1));
      assign(fh1, src_src1);
      // the edge and soucre are in the same face, so we set halfedge-target
      // information to true
      res_h->twin()->set_is_equal_aux_data_in_target(0, true);
      res_h->twin()->set_has_equal_aux_data_in_target(0, !fh1->has_no_data());

//      set_halfedge_target_info(res_h->twin(), 0, true);
    }
  }

protected:

  template <class Halfedge_handle_t>
  void copy_halfedge_target_info(Halfedge_handle_t from, Res_halfedge_handle to, unsigned int id)
  {
    to->set_is_equal_aux_data_in_target(id, from->get_is_equal_data_in_target());
    to->set_has_equal_aux_data_in_target(id, from->get_has_equal_data_in_target());    
  }
  void set_halfedge_target_info(Res_halfedge_handle to, unsigned int id, bool info)
  {
    to->set_is_equal_aux_data_in_target(id, info);
    to->set_has_equal_aux_data_in_target(id, info);
  }

  template <class Halfedge_handle_t>
  void copy_halfedge_target_info_from_halfedge_face_info(Halfedge_handle_t from,
                                                         Res_halfedge_handle to,
                                                         unsigned int id)
  {
    to->set_is_equal_aux_data_in_target(id, from->get_is_equal_data_in_face());
    to->set_has_equal_aux_data_in_target(id, from->get_has_equal_data_in_face());
  }
  template <class Vertex_handle_t>
  void copy_halfedge_target_info_from_vertex_face_info(Vertex_handle_t from,
                                                       Res_halfedge_handle to,
                                                       unsigned int id)
  {
    to->set_is_equal_aux_data_in_target(id, from->get_is_equal_data_in_face());
    to->set_has_equal_aux_data_in_target(id, from->get_has_equal_data_in_face());
  }
  
  Minimization_diagram_2& m_1;
  Minimization_diagram_2& m_2;
  Minimization_diagram_2& m_result;  
};

CGAL_END_NAMESPACE

#endif


