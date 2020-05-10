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

#ifndef CGAL_DRAW_POLYGON_WITH_HOLES_2_H
#define CGAL_DRAW_POLYGON_WITH_HOLES_2_H

#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
\ingroup PkgDrawPolygonWithHoles2

opens a new window and draws `aph`, an instance of the `CGAL::Polygon_with_holes_2` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam PH an instance of the `CGAL::Polygon_with_holes_2` class.
\param aph the polygon with holes to draw.

*/
template<class PH>
void draw(const PH& aph);

} /* namespace CGAL */
#endif

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Random.h>

namespace CGAL
{

// Viewer class for Polygon_2
template<class P2>
class SimplePolygonWithHoles2ViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt      Base;
  typedef typename P2::General_polygon_2::Point_2 Point;

public:
  /// Construct the viewer.
  /// @param ap2 the polygon to view
  /// @param title the title of the window
  SimplePolygonWithHoles2ViewerQt(QWidget* parent, const P2& ap2,
                                  const char* title="Basic Polygon_with_holes_2 Viewer") :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false),
    p2(ap2)
  {
    compute_elements();
  }

protected:
  void compute_one_loop_elements(const typename P2::General_polygon_2& p, bool hole)
  {
    if (hole)
    { add_point_in_face(p.vertex(p.size()-1)); }

    typename P2::General_polygon_2::Vertex_const_iterator prev;
    for (typename P2::General_polygon_2::Vertex_const_iterator i=p.vertices_begin();
         i!=p.vertices_end(); ++i)
    {
      add_point(*i);         // Add vertex
      if (i!=p.vertices_begin())
      { add_segment(*prev, *i); } // Add segment with previous point
      add_point_in_face(*i); // Add point in face
      prev=i;
    }

    // Add the last segment between the last point and the first one
    add_segment(*prev, *(p.vertices_begin()));
  }

  void compute_elements()
  {
    clear();

    if (p2.outer_boundary().is_empty()) return;

    CGAL::Color c(75,160,255);
    face_begin(c);

    compute_one_loop_elements(p2.outer_boundary(), false);

    for (typename P2::Hole_const_iterator it=p2.holes_begin(); it!=p2.holes_end(); ++it)
    {
      compute_one_loop_elements(*it, true);
      add_point_in_face(p2.outer_boundary().vertex(p2.outer_boundary().size()-1));
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
void draw(const CGAL::Polygon_with_holes_2<T, C>& ap2,
          const char* title="Polygon_with_holes_2 Basic Viewer")
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
    SimplePolygonWithHoles2ViewerQt<CGAL::Polygon_with_holes_2<T, C> >
      mainwindow(app.activeWindow(), ap2, title);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYGON_WITH_HOLES_2_H
