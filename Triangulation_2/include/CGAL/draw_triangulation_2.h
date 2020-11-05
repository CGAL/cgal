// Copyright(c) 2018  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_T2_H
#define CGAL_DRAW_T2_H

#include <CGAL/license/Triangulation_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Triangulation_2.h>
#include <CGAL/Random.h>

namespace CGAL
{

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorT2
{
  template<typename T2>
  static CGAL::Color run(const T2&,
                         const typename T2::Finite_faces_iterator fh)
  {
    CGAL::Random random((unsigned int)(std::size_t)(&*fh));
    return get_random_color(random);
  }
};

// Viewer class for T2
template<class T2, class ColorFunctor>
class SimpleTriangulation2ViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt                    Base;
  typedef typename T2::Vertex_handle         Vertex_const_handle;
  typedef typename T2::Finite_edges_iterator Edge_const_handle;
  typedef typename T2::Finite_faces_iterator Facet_const_handle;
  typedef typename T2::Point                 Point;

public:
  /// Construct the viewer.
  /// @param at2 the t2 to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimpleTriangulation2ViewerQt(QWidget* parent, const T2& at2,
                               const char* title="Basic T2 Viewer",
                               bool anofaces=false,
                               const ColorFunctor& fcolor=ColorFunctor()) :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false),
    t2(at2),
    m_nofaces(anofaces),
    m_fcolor(fcolor)
  {
    compute_elements();
  }

protected:
  void compute_face(Facet_const_handle fh)
  {
    CGAL::Color c=m_fcolor.run(t2, fh);
    face_begin(c);

    add_point_in_face(fh->vertex(0)->point());
    add_point_in_face(fh->vertex(1)->point());
    add_point_in_face(fh->vertex(2)->point());

    face_end();
  }

  void compute_edge(Edge_const_handle eh)
  {
    add_segment(eh->first->vertex(eh->first->cw(eh->second))->point(),
                eh->first->vertex(eh->first->ccw(eh->second))->point());
  }

  void compute_vertex(Vertex_const_handle vh)
  { add_point(vh->point()); }

  void compute_elements()
  {
    clear();

    if (!m_nofaces)
    {
      for (typename T2::Finite_faces_iterator it=t2.finite_faces_begin();
           it!=t2.finite_faces_end(); ++it)
      { compute_face(it); }
    }

    for (typename T2::Finite_edges_iterator it=t2.finite_edges_begin();
         it!=t2.finite_edges_end(); ++it)
    { compute_edge(it); }

    for (typename T2::Finite_vertices_iterator it=t2.finite_vertices_begin();
         it!=t2.finite_vertices_end(); ++it)
    { compute_vertex(it); }
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
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
  const T2& t2;
  bool m_nofaces;
  const ColorFunctor& m_fcolor;
};

// Specialization of draw function.
#define CGAL_T2_TYPE CGAL::Triangulation_2<Gt, Tds>

template<class Gt, class Tds>
void draw(const CGAL_T2_TYPE& at2,
          const char* title="Triangulation_2 Basic Viewer",
          bool nofill=false)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    int argc=1;
    const char* argv[2]={"t2_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    DefaultColorFunctorT2 fcolor;
    SimpleTriangulation2ViewerQt<CGAL_T2_TYPE, DefaultColorFunctorT2>
      mainwindow(app.activeWindow(), at2, title, nofill, fcolor);
    mainwindow.show();
    app.exec();
  }
}

#undef CGAL_T2_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_T2_H
