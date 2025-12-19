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
//             Baruch Zukerman        <baruchzu@post.tau.ac.il>
//             Efi Fogel              <efif@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_OVERLAY_2_H
#define CGAL_ENVELOPE_OVERLAY_2_H

#include <CGAL/license/Envelope_3.h>


#include <iostream>

#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Envelope_3/Envelope_overlay_functor.h>

namespace CGAL {

template <typename MinimizationDiagram_2,
          typename OverlayFunctor =
            Envelope_overlay_functor<MinimizationDiagram_2>>
class Envelope_overlay_2 {
public:
  using Minimization_diagram_2 = MinimizationDiagram_2;

  using Face_handle = typename Minimization_diagram_2::Face_handle;
  using Face_iterator = typename Minimization_diagram_2::Face_iterator;

  using Vertex_handle = typename Minimization_diagram_2::Vertex_handle;
  using Vertex_iterator = typename Minimization_diagram_2::Vertex_iterator;

  using Halfedge_handle = typename Minimization_diagram_2::Halfedge_handle;
  using Halfedge_iterator = typename Minimization_diagram_2::Halfedge_iterator;

  using Overlay_functor = OverlayFunctor;

protected:
  using Traits = typename Minimization_diagram_2::Geometry_traits_2;
  using Xy_monotone_surface_3 = typename Traits::Xy_monotone_surface_3;

public:
  void operator()(Minimization_diagram_2& md1, Minimization_diagram_2& md2,
                  Minimization_diagram_2& result) {
    CGAL_assertion(md1.is_valid());
    CGAL_assertion(md2.is_valid());

    Overlay_functor overlay_func(md1, md2, result);
    overlay(md1, md2, result, overlay_func);

    CGAL_assertion_code(post_test_assertions(result));
  }

public:
  /* void print_face(Face_handle fh) {
   *   std::cout << (fh->is_unbounded() ? "unbounded" : "bounded");
   *
   *   if (fh->is_env_set()) {
   *     std::cout << " #data= " << fh->env_data_size();
   *     if (fh->env_data_size() > 0)
   *       std::cout << " data= " << fh->env_data_front();
   *   }
   *
   *   if (fh->aux_is_set(0)) {
   *     std::cout << " #data1= " << number_of_aux_data_objects(fh, 0);
   *     if (number_of_aux_data_objects(fh, 0)>0)
   *       std::cout << " data#1= " << aux_data(fh, 0);
   *   }
   *   if (fh->aux_is_set(1)) {
   *     std::cout << " #data2= " << number_of_aux_data_objects(fh, 1);
   *     if (number_of_aux_data_objects(fh, 1)>0)
   *       std::cout << " data#2= " << aux_data(fh, 1);
   *   }
   *   std::cout << std::endl;
   * }
   *
   * // print the aux data in the faces of md
   * void print_faces(Minimization_diagram_2& md) {
   *   Face_iterator fit = md.faces_begin();
   *   for(; fit != md.faces_end(); ++fit) {
   *     Face_handle fh = fit;
   *     print_face(fh);
   *   }
   *   std::cout << std::endl;
   * }
   *
   * void print_vertices(Minimization_diagram_2& md) {
   *   Vertex_iterator vit = md.vertices_begin();
   *   for(; vit != md.vertices_end(); ++vit) {
   *     Vertex_handle vh = vit;
   *     std::cout << vh->point();
   *
   *      if (vh->is_env_set()) {
   *       std::cout << " #data= " << vh->env_data_size();
   *       if (vh->env_data_size() > 0)
   *         std::cout << " data= " << vh->env_data_front();
   *     }
   *
   *     if (vh->aux_is_set(0)) {
   *       std::cout << " #data1= " << number_of_aux_data_objects(vh, 0);
   *       if (number_of_aux_data_objects(vh, 0)>0)
   *         std::cout << " data#1= " << aux_data(vh, 0);
   *     }
   *     if (vh->aux_is_set(1)) {
   *       std::cout << " #data2= " << number_of_aux_data_objects(vh, 1);
   *       if (number_of_aux_data_objects(vh, 1)>0)
   *         std::cout << " data#2= " << aux_data(vh, 1);
   *     }
   *     std::cout << std::endl;
   *   }
   *   std::cout << std::endl;
   * }
   *
   * void print_edges(Minimization_diagram_2& md) {
   *   Halfedge_iterator hit = md.halfedges_begin();
   *   for(; hit != md.halfedges_end(); ++hit, ++hit) {
   *     Halfedge_handle hh = hit;
   *     std::cout << hh->curve();
   *
   *     if (hh->is_env_set()) {
   *       std::cout << " #data= " << hh->env_data_size();
   *       if (hh->env_data_size() > 0)
   *         std::cout << " data= " << hh->env_data_front();
   *     }
   *
   *     if (hh->aux_is_set(0)) {
   *       std::cout << " #data1= " << number_of_aux_data_objects(hh, 0);
   *       if (number_of_aux_data_objects(hh, 0)>0)
   *         std::cout << " data#1= " << aux_data(hh, 0);
   *     }
   *     if (hh->aux_is_set(1)) {
   *       std::cout << " #data2= " << number_of_aux_data_objects(hh, 1);
   *
   *       if (number_of_aux_data_objects(hh, 1)>0)
   *         std::cout << " data#2= " << aux_data(hh, 1);
   *     }
   *     std::cout << std::endl;
   *   }
   *   std::cout << std::endl;
   * }
   */

  void post_test_assertions(Minimization_diagram_2& md) {
    // check that all data is filled in result
    for (auto fi = md.faces_begin(); fi != md.faces_end(); ++fi) {
      Face_handle fh = fi;
      CGAL_assertion_msg(fh->aux_is_set(0),
                         "data from md1 on face is not set");
      CGAL_assertion_msg(fh->aux_is_set(1),
                         "data from md2 on face is not set");
    }

    for (auto hi = md.halfedges_begin(); hi != md.halfedges_end(); ++hi) {
      Halfedge_handle hh = hi;
      CGAL_assertion_msg(hh->aux_is_set(0),
                         "data from md1 on halfedge is not set");
      CGAL_assertion_msg(hh->aux_is_set(1),
                         "data from md2 on halfedge is not set");
    }

    for (auto vi = md.vertices_begin(); vi != md.vertices_end(); ++vi) {
      Vertex_handle vh = vi;
      CGAL_assertion_msg(vh->aux_is_set(0),
                         "data from md1 on vertex is not set");
      CGAL_assertion_msg(vh->aux_is_set(1),
                         "data from md2 on vertex is not set");
    }
  }

protected:
  // helper methods
  template <typename FeatureHandle>
  Xy_monotone_surface_3 aux_data(FeatureHandle fh, unsigned int id) {
    const Object& o = fh->aux_source(id);
    Xy_monotone_surface_3 data;

    Halfedge_handle h;
    Vertex_handle v;
    Face_handle f;
    if (assign(v, o)) data = v->env_data_front();
    else if (assign(h, o)) data = h->env_data_front();
    else {
      CGAL_assertion(assign(f, o));
      assign(f, o);
      data = f->env_data_front();
    }
    return data;
  }

  template <typename FeatureHandle>
  int number_of_aux_data_objects(FeatureHandle fh, unsigned int id) {
    const Object& o = fh->aux_source(id);
    int data;

    Halfedge_handle h;
    Vertex_handle v;
    Face_handle f;
    if (assign(v, o)) data = v->env_data_size();
    else if (assign(h, o)) data = h->env_data_size();
    else {
      CGAL_assertion(assign(f, o));
      assign(f, o);
      data = f->env_data_size();
    }
    return data;
  }

};

} //namespace CGAL

#endif
