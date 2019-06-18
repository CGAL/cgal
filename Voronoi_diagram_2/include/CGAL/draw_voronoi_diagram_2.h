// Copyright(c) 2018  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Jasmeet Singh <jasmeet.singh.mec11@iitbhu.ac.in>

#ifndef CGAL_DRAW_VORONOI_H
#define CGAL_DRAW_VORONOI_H

#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/license/Voronoi_diagram_2.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/Converter.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Voronoi_diagram_2/Face.h>
#include <CGAL/Voronoi_diagram_2/Handle_adaptor.h>
#include <CGAL/Voronoi_diagram_2/Vertex.h>
#include <CGAL/Voronoi_diagram_2/basic.h>

namespace CGAL {

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorV2 {
  template <typename V2>
  static CGAL::Color run(const V2 &, const typename V2::Face_iterator fh) {
    CGAL::Random random((unsigned int)(std::size_t)(&*fh));
    return get_random_color(random);
  }
};

// Viewer for Voronoi diagram
template <class V2, class ColorFunctor>
class SimpleVoronoiDiagram2ViewerQt : public Basic_viewer_qt {
  typedef Basic_viewer_qt Base;
  typedef typename V2::Halfedge_iterator Halfedge_const_handle;
  typedef typename V2::Face_iterator Face_const_handle;
  typedef typename V2::Vertex_iterator Vertex_const_handle;
  typedef typename V2::Point_2 Point;
  typedef typename V2::Delaunay_vertex_handle Delaunay_vertex_const_handle;
  typedef typename V2::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  typedef typename V2::Delaunay_geom_traits Delaunay_geom_traits;

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Face<V2> Face;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Handle_adaptor<Face> Face_handle;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Vertex<V2> Vertex;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Handle_adaptor<Vertex> Vertex_handle;
  typedef Triangulation_cw_ccw_2 CW_CCW_2;

public:
  /// Construct the viewer.
  /// @param av2 the voronoi diagram to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this
  /// can be
  ///        useful for very big object where this time could be long)
  SimpleVoronoiDiagram2ViewerQt(QWidget *parent, const V2 &av2,
                                const char *title = "Basic Voronoi Viewer",
                                bool anofaces = false,
                                const ColorFunctor &fcolor = ColorFunctor())
      : // First draw: vertices; half-edges; faces; multi-color; no inverse
        // normal
        Base(parent, title, true, true, true, false, false), v2(av2),
        m_nofaces(anofaces), m_fcolor(fcolor) {
    compute_elements();
  }

protected:
  void compute_face(Face_const_handle fh) {
    //    CGAL::Color c=m_fcolor.run(v2, fh);
    //    face_begin(c);
    //    Ccb_halfedge_circulator ec_start = fh->ccb();
    //    Ccb_halfedge_circulator ec = ec_start;

    //    do{
    //        if( ec->has_source() )
    //            add_point_in_face(ec->source()->point());
    //        else if(ec->has_target())
    //            add_point_in_face(ec->target()->point());
    //    } while(++ec != ec_start);

    //    face_end();
  }

  void compute_ray_points(Halfedge_const_handle he) {
    if (he->is_segment()) {
      add_segment(he->source()->point(), he->target()->point());
    } else if (he->is_ray()) {
      Delaunay_vertex_const_handle v1 = he->up();
      Delaunay_vertex_const_handle v2 = he->down();

      Kernel::Vector_2 direction(v1->point().y() - v2->point().y(),
                                 v2->point().x() - v1->point().x());
      Kernel::Point_2 end_point;

      if (he->has_source()) {
        end_point = he->source()->point();
        add_ray_points(end_point, direction);
      }
    }
  }

  void compute_rays(Halfedge_const_handle he) {
    if (he->is_ray()) {
      Delaunay_vertex_const_handle v1 = he->up();
      Delaunay_vertex_const_handle v2 = he->down();

      Kernel::Vector_2 direction(v1->point().y() - v2->point().y(),
                                 v2->point().x() - v1->point().x());
      Kernel::Point_2 end_point;
      if (he->has_source()) {
        end_point = he->source()->point();
        std::cout << "Bounding_box" << m_bounding_box << std::endl;
        add_ray_segment(end_point, direction);
      }
    }
  }

  void compute_vertex(Vertex_const_handle vh) { add_point(vh->point()); }

  void compute_elements() {
    clear();

    if (!m_nofaces) {
      for (typename V2::Face_iterator it = v2.faces_begin();
           it != v2.faces_end(); ++it) {
        compute_face(it);
      }
    }
    //    for(Delaunay_vertex_const_handle it =
    //    v2.dual().finite_vertices_begin();
    //        it!=v2.dual().finite_vertices_end(); ++it)
    //    { compute_vertex(it);}

    for (typename V2::Halfedge_iterator it = v2.halfedges_begin();
         it != v2.halfedges_end(); ++it) {
      compute_ray_points(it);
    }

    for (typename V2::Halfedge_iterator it = v2.halfedges_begin();
         it != v2.halfedges_end(); ++it) {
      compute_rays(it);
    }

    for (typename V2::Vertex_iterator it = v2.vertices_begin();
         it != v2.vertices_end(); ++it) {
      compute_vertex(it);
    }
  }

  virtual void keyPressEvent(QKeyEvent *e) {
    // Test key pressed:
    //    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    //    if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton)) { ... }

    // Call: * compute_elements() if the model changed, followed by
    //       * redraw() if some viewing parameters changed that implies some
    //                  modifications of the buffers
    //                  (eg. type of normal, color/mono)
    //       * update() just to update the drawing

    // Call the base method to process others/classicals key
    Base::keyPressEvent(e);
  }

protected:
  CGAL::Qt::Converter<Delaunay_geom_traits> convert;
  const V2 &v2;
  bool m_nofaces;
  const ColorFunctor &m_fcolor;
};

template <class V2, class ColorFunctor>
void draw(const V2 &av2, const char *title, bool nofill,
          const ColorFunctor &fcolor) {
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite = true;
#else
  bool cgal_test_suite = false;
#endif

  if (!cgal_test_suite) {
    int argc = 1;
    const char *argv[2] = {"v2_viewer", "\0"};
    QApplication app(argc, const_cast<char **>(argv));
    SimpleVoronoiDiagram2ViewerQt<V2, ColorFunctor> mainwindow(
        app.activeWindow(), av2, title, nofill, fcolor);
    mainwindow.show();
    app.exec();
  }
}

template <class V2> void draw(const V2 &av2, const char *title, bool nofill) {
  DefaultColorFunctorV2 c;
  draw(av2, title, nofill, c);
}

template <class V2> void draw(const V2 &av2, const char *title) {
  draw(av2, title, false);
}

template <class V2> void draw(const V2 &av2) {
  draw(av2, "Basic Voronoi Diagram Viewer");
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_VORONOI_H
