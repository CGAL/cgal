// Copyright (c) 2016  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_POINT_SET_3_H
#define CGAL_DRAW_POINT_SET_3_H

#include <CGAL/license/Point_set_3.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
\ingroup PkgDrawPointSet3D

opens a new window and draws `aps`, an instance of the `CGAL::Point_set_3` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam PS an instance of the `CGAL::Point_set_3` class.
\param aps the point set to draw.

*/
template<class PS>
void draw(const PS& aps);

} /* namespace CGAL */
#endif

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Point_set_3.h>
#include <CGAL/Random.h>

namespace CGAL
{

// Viewer class for Point_set
template<class PointSet>
class SimplePointSetViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt Base;
  typedef typename PointSet::Point_map::value_type Point;

public:
  /// Construct the viewer.
  /// @param apointset the point set to view
  /// @param title the title of the window
  SimplePointSetViewerQt(QWidget* parent,
                         const PointSet& apointset, const char* title="") :
    // First draw: vertices; no-edge, no-face; mono-color; no inverse normal
    Base(parent, title, true, false, false, true, false),
    pointset(apointset)
  {
    compute_elements();
  }

protected:
  void compute_vertex(const Point& p)
  {
    add_point(p);
    // We can use add_point(p, c) with c a CGAL::Color to add a colored point
  }

  void compute_elements()
  {
    clear();

    for (typename PointSet::const_iterator it=pointset.begin();
          it!=pointset.end(); ++it)
    { compute_vertex(pointset.point(*it)); }
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    // const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    Base::keyPressEvent(e);
  }

protected:
  const PointSet& pointset;
};

// Specialization of draw function.
template<class P, class V>
void draw(const Point_set_3<P, V>& apointset,
          const char* title="Point_set_3 Basic Viewer")
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    int argc=1;
    const char* argv[2]={"point_set_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    SimplePointSetViewerQt<Point_set_3<P, V> > mainwindow(app.activeWindow(),
                                                          apointset,
                                                          title);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POINT_SET_3_H
