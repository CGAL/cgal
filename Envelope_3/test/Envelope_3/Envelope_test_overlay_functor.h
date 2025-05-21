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
// Author(s) : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//             Efi Fogel              <efif@post.tau.ac.il>

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
template <typename MinimizationDiagram_2>
class Envelope_test_overlay_functor {
public:
  using Minimization_diagram_2 = MinimizationDiagram_2;

private:
  using Md2 = Minimization_diagram_2;

public:
  using Face_handle1 = typename Md2::Face_const_handle;
  using Face_handle2 = typename Md2::Face_const_handle;

  using Vertex_handle1 = typename Md2::Vertex_const_handle;
  using Vertex_handle2 = typename Md2::Vertex_const_handle;

  using Halfedge_handle1 = typename Md2::Halfedge_const_handle;
  using Halfedge_handle2 = typename Md2::Halfedge_const_handle;

  using Res_face_handle = typename Md2::Face_handle;
  using Res_halfedge_handle = typename Md2::Halfedge_handle;
  using Res_vertex_handle = typename Md2::Vertex_handle;

  Envelope_test_overlay_functor(Md2&, Md2&, Md2&) {}

  void create_face (Face_handle1 f1, Face_handle2 f2, Res_face_handle res_f) {
    res_f->set_aux_source(0, f1);
    res_f->set_aux_source(1, f2);
    assert_msg(f1->is_equal_env_data(f2->begin_env_data(), f2->end_env_data()),
               "data different over face");
  }

  void create_vertex(Halfedge_handle1 h1, Halfedge_handle2 h2,
                     Res_vertex_handle res_v) {
    res_v->set_aux_source(0, h1);
    res_v->set_aux_source(1, h2);
    assert_msg(h1->is_equal_env_data(h2->begin_env_data(), h2->end_env_data()),
               "data different over vertex");
  }

  void create_vertex(Vertex_handle1 v1, Vertex_handle2 v2,
                     Res_vertex_handle res_v) {
    res_v->set_aux_source(0, v1);
    res_v->set_aux_source(1, v2);
    assert_msg(v1->is_equal_env_data(v2->begin_env_data(), v2->end_env_data()),
               "data different over vertex");
  }

  void create_vertex(Vertex_handle1 v1, Halfedge_handle2 h2,
                     Res_vertex_handle res_v) {
    res_v->set_aux_source(0, v1);
    res_v->set_aux_source(1, h2);
    assert_msg(v1->is_equal_env_data(h2->begin_env_data(), h2->end_env_data()),
               "data different over vertex");
  }

  void create_vertex(Halfedge_handle1 h1, Vertex_handle2 v2,
                     Res_vertex_handle res_v) {
    res_v->set_aux_source(0, h1);
    res_v->set_aux_source(1, v2);
    assert_msg(h1->is_equal_env_data(v2->begin_env_data(), v2->end_env_data()),
               "data different over vertex");
  }

  void create_vertex(Face_handle1 f1, Vertex_handle2 v2,
                     Res_vertex_handle res_v) {
    res_v->set_aux_source(0, f1);
    res_v->set_aux_source(1, v2);
    assert_msg(f1->is_equal_env_data(v2->begin_env_data(), v2->end_env_data()),
               "data different over vertex");
  }

  void create_vertex(Vertex_handle1 v1, Face_handle2 f2,
                     Res_vertex_handle res_v) {
    res_v->set_aux_source(0, v1);
    res_v->set_aux_source(1, f2);
    assert_msg(v1->is_equal_env_data(f2->begin_env_data(), f2->end_env_data()),
               "data different over vertex");
  }

  void create_edge(Halfedge_handle1 h1, Halfedge_handle2 h2,
                   Res_halfedge_handle res_h) {
    res_h->set_aux_source(0, h1);
    res_h->set_aux_source(1, h2);

    res_h->twin()->set_aux_source(0, h1->twin());
    res_h->twin()->set_aux_source(1, h2->twin());

    assert_msg(h1->is_equal_env_data(h2->begin_env_data(), h2->end_env_data()),
               "data different over edge");
  }

  void create_edge(Halfedge_handle1 h1, Face_handle2 f2,
                   Res_halfedge_handle res_h) {
    res_h->set_aux_source(0, h1);
    res_h->set_aux_source(1, f2);

    res_h->twin()->set_aux_source(0, h1->twin());
    res_h->twin()->set_aux_source(1, f2);

    assert_msg(h1->is_equal_env_data(f2->begin_env_data(), f2->end_env_data()),
               "data different over edge");
  }

  void create_edge(Face_handle1 f1, Halfedge_handle2 h2,
                   Res_halfedge_handle res_h) {
    res_h->set_aux_source(0, f1);
    res_h->set_aux_source(1, h2);

    res_h->twin()->set_aux_source(0, f1);
    res_h->twin()->set_aux_source(1, h2->twin());
    assert_msg(f1->is_equal_env_data(h2->begin_env_data(), h2->end_env_data()),
               "data different over edge");
  }
};

} //namespace CGAL

#endif
