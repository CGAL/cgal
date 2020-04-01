// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_POLYGON_2_H
#define CGAL_DRAW_POLYGON_2_H

#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
\ingroup PkgDrawPolygon2

opens a new window and draws `ap`, an instance of the `CGAL::Polygon_2` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam P an instance of the `CGAL::Polygon_2` class.
\param ap the polygon to draw.

*/
template<class P>
void draw(const P& ap);

} /* namespace CGAL */

#endif

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Polygon_2.h>
#include <CGAL/Random.h>

namespace CGAL
{

// Viewer class for Polygon_2
template<class P2>
class SimplePolygon2ViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt      Base;
  typedef typename P2::Point_2 Point;

public:
  /// Construct the viewer.
  /// @param ap2 the polygon to view
  /// @param title the title of the window
  SimplePolygon2ViewerQt(QWidget* parent, const P2& ap2,
                         const char* title="Basic Polygon_2 Viewer") :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false),
    p2(ap2)
  {
    compute_elements();
  }

protected:
  void compute_elements()
  {
    clear();

    if (p2.is_empty()) return;

    Point prev=p2.vertex(p2.size()-1);

    CGAL::Color c(75,160,255);
    face_begin(c);

    for (typename P2::Vertex_const_iterator i=p2.vertices_begin();
         i!=p2.vertices_end(); ++i)
    {
      add_point(*i);         // Add vertex
      add_segment(prev, *i); // Add segment with previous point
      add_point_in_face(*i); // Add point in face
      prev=*i;
    }

    face_end();
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
  const P2& p2;
};

// Specialization of draw function.
template<class T, class C>
void draw(const CGAL::Polygon_2<T, C>& ap2,
          const char* title="Polygon_2 Basic Viewer")
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
    SimplePolygon2ViewerQt<CGAL::Polygon_2<T, C> >
      mainwindow(app.activeWindow(), ap2, title);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYGON_2_H
