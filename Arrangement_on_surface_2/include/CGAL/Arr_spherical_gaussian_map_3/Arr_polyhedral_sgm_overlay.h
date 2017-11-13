// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University(Israel).
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
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_H
#define CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_H

#include <CGAL/license/Arrangement_on_surface_2.h>


// #define CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG 1

#include <iostream>
#include <CGAL/use.h>

namespace CGAL {

template <class Sgm>
class Arr_polyhedral_sgm_overlay {
private:
  typedef typename Sgm::Point_3                         Point_3;
  typedef typename Sgm::Vector_3                        Vector_3;
  
public:
  typedef typename Sgm::Face_handle                     Face_handle;
  typedef typename Sgm::Vertex_handle                   Vertex_handle;
  typedef typename Sgm::Halfedge_handle                 Halfedge_handle;

  typedef typename Sgm::Face_const_handle               Face_const_handle;
  typedef typename Sgm::Vertex_const_handle             Vertex_const_handle;
  typedef typename Sgm::Halfedge_const_handle           Halfedge_const_handle;

  typedef typename Sgm::Ccb_halfedge_const_circulator 
    Ccb_halfedge_const_circulator;

  typedef typename Sgm::Halfedge_around_vertex_const_circulator
    Arr_halfedge_around_vertex_const_circulator;
  
  /*! 1 */
  void create_face(Face_const_handle f1, Face_const_handle f2, Face_handle f)
  {
    const Point_3 & p1 = f1->point();
    const Point_3 & p2 = f2->point();
    Vector_3 v1(ORIGIN, p1);
    Point_3 p = p2 + v1;
    f->set_point(p);
#ifdef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG
    std::cout << "create_face(f1, f2) "
              << p
              << std::endl;
    std::cout << "  Outer CCB:" << std::endl;
    typename Sgm::Outer_ccb_iterator oit;
    for (oit = f->outer_ccbs_begin(); oit != f->outer_ccbs_end(); ++oit) {
      typename Sgm::Halfedge_iterator first = *oit;
      typename Sgm::Halfedge_iterator curr = first;
      do {
        std::cout << "  " << curr->curve() << std::endl;
        curr = curr->next();
      } while (curr != first);
    }

    std::cout << "  Inner CCB:" << std::endl;
    typename Sgm::Inner_ccb_iterator iit;
    for (iit = f->inner_ccbs_begin(); iit != f->inner_ccbs_end(); ++iit) {
      typename Sgm::Halfedge_iterator first = *iit;
      typename Sgm::Halfedge_iterator curr = first;
      do {
        std::cout << "  " << curr->curve() << std::endl;
        curr = curr->next();
      } while (curr != first);
    }
#endif
  }

  /*! 2 */
  void create_vertex(Halfedge_const_handle /*h1*/, Halfedge_const_handle /*h2*/,
                     Vertex_handle v)
  {
#ifdef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG
    std::cout << "create_vertex(h1, h2)"
              << " " << v->point()
              << std::endl;
#else
    CGAL_USE(v);
#endif
  }

  /*! 3 */
  void create_vertex(Vertex_const_handle /*v1*/, Vertex_const_handle /*v2*/,
                     Vertex_handle v)
  {
#ifdef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG
    std::cout << "create_vertex(v1, v2)"
              << " " << v->point()
              << std::endl;
#else
    CGAL_USE(v);
#endif
  }

  /*! 4 */
  void create_vertex(Vertex_const_handle /*v1*/, Halfedge_const_handle /*h2*/,
                     Vertex_handle v)
  {
#ifdef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG
    std::cout << "create_vertex(v1, h2)"
              << " " << v->point()
              << std::endl;
#else
    CGAL_USE(v);
#endif
  }

  /*! 5 */
  void create_vertex(Halfedge_const_handle /*h1*/, Vertex_const_handle /*v2*/,
                     Vertex_handle v)
  {
#ifdef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG
    std::cout << "create_vertex(h1, v2)"
              << " " << v->point()
              << std::endl;
#else
    CGAL_USE(v);
#endif
  }

  /*! 6 */
  void create_vertex(Face_const_handle /*f1*/, Vertex_const_handle /*v2*/,
                     Vertex_handle v)
  {
#ifdef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG
    std::cout << "create_vertex(f1, v2)"
              << " " << v->point()
              << std::endl;
#else
    CGAL_USE(v);
#endif
  }

  /*! 7 */
  void create_vertex(Vertex_const_handle /*v1*/, Face_const_handle /*f2*/,
                     Vertex_handle v)
  {
#ifdef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG
    std::cout << "create_vertex(v1, f2)"
              << " " << v->point()
              << std::endl;
#else
    CGAL_USE(v);
#endif
  }

  /*! 8 */
  void create_edge(Halfedge_const_handle /*h1*/, Halfedge_const_handle /*h2*/,
                   Halfedge_handle h)
  {
    h->add_arr(0);
    h->add_arr(1);
    h->twin()->add_arr(0);
    h->twin()->add_arr(1);
#ifdef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG
    std::cout << "create_edge(h1, h2)"
              << " " << h->curve()
              << std::endl;
#endif
  }

  /*! 9 */
  void create_edge(Halfedge_const_handle /*h1*/, Face_const_handle /*f2*/,
                   Halfedge_handle h)
  {
    h->add_arr(0);
    h->twin()->add_arr(0);
#ifdef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG
    std::cout << "create_edge(h1, f2)"
              << " " << h->curve()
              << std::endl;
#endif
  }

  /*! 10 */
  void create_edge(Face_const_handle /*f1*/, Halfedge_const_handle /*h2*/,
                   Halfedge_handle h)
  {
    h->add_arr(1);
    h->twin()->add_arr(1);
#ifdef CGAL_ARR_POLYHEDRAL_SGM_OVERLAY_DEBUG
    std::cout << "create_edge(f1, h2)"
              << " " << h->curve()
              << std::endl;
#endif
  }
};

} //namespace CGAL

#endif
