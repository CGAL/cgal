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
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef ENVELOPE_TEST_OVERLAY_FUNCTOR_H
#define ENVELOPE_TEST_OVERLAY_FUNCTOR_H

#ifndef assert_msg
#ifndef NDEBUG
#  define assert_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#else
#  define assert_msg(EX,MSG) (static_cast<void>(0))
#endif
#endif

#include <iostream>

namespace CGAL {

// this overlay functor compares the data over the 2 features that create new
// features in the new map
template <class MinimizationDiagram_2>
class Envelope_test_overlay_functor
{
public:
  typedef MinimizationDiagram_2                                  Minimization_diagram_2;
  
  typedef typename Minimization_diagram_2::Face_const_handle     Face_handle1;
  typedef typename Minimization_diagram_2::Face_const_handle     Face_handle2;

  typedef typename Minimization_diagram_2::Vertex_const_handle   Vertex_handle1;
  typedef typename Minimization_diagram_2::Vertex_const_handle   Vertex_handle2;

  typedef typename Minimization_diagram_2::Halfedge_const_handle Halfedge_handle1;
  typedef typename Minimization_diagram_2::Halfedge_const_handle Halfedge_handle2;

  typedef typename Minimization_diagram_2::Face_handle           Res_face_handle;
  typedef typename Minimization_diagram_2::Halfedge_handle       Res_halfedge_handle;
  typedef typename Minimization_diagram_2::Vertex_handle         Res_vertex_handle;
  

  Envelope_test_overlay_functor(Minimization_diagram_2& ,
				Minimization_diagram_2& ,
				Minimization_diagram_2& )
  {}

  void create_face (Face_handle1 f1, Face_handle2 f2, Res_face_handle res_f)
  {
    res_f->set_aux_source(0, f1);
    res_f->set_aux_source(1, f2);
    assert_msg(f1->is_equal_data(f2->begin_data(), f2->end_data()),
                       "data different over face");
  }

  void create_vertex(Halfedge_handle1 h1,
                     Halfedge_handle2 h2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, h1);
    res_v->set_aux_source(1, h2);
    assert_msg(h1->is_equal_data(h2->begin_data(), h2->end_data()),
                       "data different over vertex");

  }

  void create_vertex(Vertex_handle1 v1,
                     Vertex_handle2 v2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, v1);
    res_v->set_aux_source(1, v2);
    assert_msg(v1->is_equal_data(v2->begin_data(), v2->end_data()),
                       "data different over vertex");
  }

  void create_vertex(Vertex_handle1 v1,
                     Halfedge_handle2 h2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, v1);
    res_v->set_aux_source(1, h2);
    assert_msg(v1->is_equal_data(h2->begin_data(), h2->end_data()),
                       "data different over vertex");
  }

  void create_vertex(Halfedge_handle1 h1,
                     Vertex_handle2 v2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, h1);
    res_v->set_aux_source(1, v2);
    assert_msg(h1->is_equal_data(v2->begin_data(), v2->end_data()),
                       "data different over vertex");
  }

  void create_vertex(Face_handle1 f1,
                     Vertex_handle2 v2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, f1);
    res_v->set_aux_source(1, v2);
    assert_msg(f1->is_equal_data(v2->begin_data(), v2->end_data()),
                       "data different over vertex");
  }

  void create_vertex(Vertex_handle1 v1,
                     Face_handle2 f2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, v1);
    res_v->set_aux_source(1, f2);
    assert_msg(v1->is_equal_data(f2->begin_data(), f2->end_data()),
                       "data different over vertex");
  }

  void create_edge(Halfedge_handle1 h1,
                   Halfedge_handle2 h2,
                   Res_halfedge_handle res_h)
  {
    res_h->set_aux_source(0, h1);
    res_h->set_aux_source(1, h2);

    res_h->twin()->set_aux_source(0, h1->twin());
    res_h->twin()->set_aux_source(1, h2->twin());

    assert_msg(h1->is_equal_data(h2->begin_data(), h2->end_data()),
                       "data different over edge");
  }

  void create_edge(Halfedge_handle1 h1,
                   Face_handle2 f2,
                   Res_halfedge_handle res_h)
  {
    res_h->set_aux_source(0, h1);
    res_h->set_aux_source(1, f2);

    res_h->twin()->set_aux_source(0, h1->twin());
    res_h->twin()->set_aux_source(1, f2);

    assert_msg(h1->is_equal_data(f2->begin_data(), f2->end_data()),
                       "data different over edge");
  }

  void create_edge(Face_handle1 f1,
                   Halfedge_handle2 h2,
                   Res_halfedge_handle res_h)
  {
    res_h->set_aux_source(0, f1);
    res_h->set_aux_source(1, h2);

    res_h->twin()->set_aux_source(0, f1);
    res_h->twin()->set_aux_source(1, h2->twin());
    assert_msg(f1->is_equal_data(h2->begin_data(), h2->end_data()),
                       "data different over edge");

  }
  
};

} //namespace CGAL

#endif


