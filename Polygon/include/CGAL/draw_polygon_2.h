// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_POLYGON_2_H
#define CGAL_DRAW_POLYGON_2_H

#include <CGAL/Qt/Basic_viewer_qt.h>

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
          const char* title="Polygon_2 Basic Viewer",
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
    SimplePolygon2ViewerQt<CGAL::Polygon_2<T, C> >
      mainwindow(app.activeWindow(), ap2, title);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYGON_2_H
