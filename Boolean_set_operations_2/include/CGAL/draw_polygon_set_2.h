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

#ifndef CGAL_DRAW_POLYGON_SET_2_H
#define CGAL_DRAW_POLYGON_SET_2_H

#include <CGAL/draw_polygon_with_holes_2.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
\ingroup PkgDrawPolygonSet2

opens a new window and draws `aps`, an instance of the `CGAL::Polygon_set_2` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam PS an instance of the `CGAL::Polygon_set_2` class.
\param aps the polygon set to draw.

*/
template<class PS>
void draw(const PS& aps);

} /* namespace CGAL */
#endif

#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/Polygon_set_2.h>

namespace CGAL
{

template<class PS2>
class SimplePolygonSet2ViewerQt :
  public SimplePolygonWithHoles2ViewerQt<typename PS2::Polygon_with_holes_2>
{
  typedef SimplePolygonWithHoles2ViewerQt<typename PS2::Polygon_with_holes_2> Base;

public:
  SimplePolygonSet2ViewerQt(QWidget* parent, const PS2& aps2,
                            const char* title="Basic Polygon_set_2 Viewer") :
    Base(parent, title)
  {
    std::vector<typename PS2::Polygon_with_holes_2> polygons;
    aps2.polygons_with_holes(std::back_inserter(polygons));

    for (typename PS2::Polygon_with_holes_2& P: polygons) {
      Base::compute_elements(P);
    }
  }
};

// Specialization of draw function.
template<class T, class C>
void draw(const CGAL::Polygon_set_2<T, C>& aps2,
          const char* title="Polygon_set_2 Basic Viewer")
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    CGAL::Qt::init_ogl_context(4,3);
    int argc=1;
    const char* argv[2]={"t2_viewer", nullptr};
    QApplication app(argc,const_cast<char**>(argv));
    SimplePolygonSet2ViewerQt<CGAL::Polygon_set_2<T, C> >
      mainwindow(app.activeWindow(), aps2, title);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYGON_SET_2_H
