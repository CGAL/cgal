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
// Author(s)     : Ken Arroyo Ohori <k.ohori@tudelft.nl>
//                 Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_MULTIPOLYGON_WITH_HOLES_2_H
#define CGAL_DRAW_MULTIPOLYGON_WITH_HOLES_2_H

#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
 * \ingroup PkgDrawMultipolygonWithHoles2
 *
 * opens a new window and draws `aph`, an instance of the
 * `CGAL::Multipolygon_with_holes_2` class. A call to this function is blocking, that
 * is the program continues as soon as the user closes the window. This function
 * requires `CGAL_Qt6`, and is only available if the macro
 * `CGAL_USE_BASIC_VIEWER` is defined.  Linking with the cmake target
 * `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition
 * `CGAL_USE_BASIC_VIEWER`.
 * \tparam PH an instance of the `CGAL::Multipolygon_with_holes_2` class.
 * \param aph the multipolygon with holes to draw.
 */

template <typename MPH>
void draw(const MPH& aph);

} /* namespace CGAL */

#endif

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Multipolygon_with_holes_2.h>

namespace CGAL {

// Viewer class for Multipolygon_with_holes_2
template <typename Multipolygon>
class Mpwh_2_basic_viewer_qt : public Basic_viewer_qt {
  using Base = Basic_viewer_qt;
  using Mpwh = Multipolygon;
  using Pwh = typename Mpwh::Polygon_with_holes_2;
  using Pgn = typename Mpwh::Polygon_2;
  using Pnt = typename Pgn::Point_2;

public:
  /// Construct the viewer.
  /// @param parent the active window to draw
  /// @param mpwh the multipolygon to view
  /// @param title the title of the window
  Mpwh_2_basic_viewer_qt(QWidget* parent, const Mpwh& mpwh,
                         const char* title = "Basic Multipolygon_with_holes_2 Viewer") :
    Base(parent, title, true, true, true, false, false),
    m_mpwh(mpwh) {
    if (mpwh.number_of_polygons_with_holes() == 0) return;

    // mimic the computation of Camera::pixelGLRatio()
    auto bbox = bounding_box();
    CGAL::qglviewer::Vec minv(bbox.xmin(), bbox.ymin(), 0);
    CGAL::qglviewer::Vec maxv(bbox.xmax(), bbox.ymax(), 0);
    auto diameter = (maxv - minv).norm();
    m_pixel_ratio = diameter / m_height;
  }

  /*! Intercept the resizing of the window.
   */
  virtual void resizeGL(int width, int height) {
    CGAL::QGLViewer::resizeGL(width, height);
    m_width = width;
    m_height = height;
    CGAL::qglviewer::Vec p;
    auto ratio = camera()->pixelGLRatio(p);
    m_pixel_ratio = ratio;
    add_elements();
  }

  /*! Obtain the pixel ratio.
   */
  double pixel_ratio() const { return m_pixel_ratio; }

  /*! Compute the bounding box.
   */
  CGAL::Bbox_2 bounding_box() {
    Bbox_2 bbox;
    if (m_mpwh.number_of_polygons_with_holes() > 0) bbox = m_mpwh.polygons_with_holes().front().outer_boundary().bbox();
    for (auto const& pwh: m_mpwh.polygons_with_holes()) {
      bbox += pwh.outer_boundary().bbox();
    }
    return bbox;
  }

  /*! Compute the elements of a multipolygon with holes.
   */
  void add_elements() {
    clear();

    for (auto const& p: m_mpwh.polygons_with_holes()) {
      CGAL::IO::Color c(rand()%255,rand()%255,rand()%255);
      face_begin(c);

      const Pnt* point_in_face;
      const auto& outer_boundary = p.outer_boundary();
      compute_loop(outer_boundary, false);
      point_in_face = &(outer_boundary.vertex(outer_boundary.size()-1));

      for (auto it = p.holes_begin(); it != p.holes_end(); ++it) {
        compute_loop(*it, true);
        add_point_in_face(*point_in_face);
      }

      face_end();
    }
  }

protected:
  /*! Compute the face
   */
  void compute_loop(const Pgn& p, bool hole) {
    if (hole) add_point_in_face(p.vertex(p.size()-1));

    auto prev = p.vertices_begin();
    auto it = prev;
    add_point(*it);
    add_point_in_face(*it);
    for (++it; it != p.vertices_end(); ++it) {
      add_segment(*prev, *it);  // add segment with previous point
      add_point(*it);
      add_point_in_face(*it);   // add point in face
      prev = it;
    }

    // Add the last segment between the last point and the first one
    add_segment(*prev, *(p.vertices_begin()));
  }

  virtual void keyPressEvent(QKeyEvent* e) {
    // Test key pressed:
    //    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    //    if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton)) { ... }

    // Call: * add_elements() if the model changed, followed by
    //       * redraw() if some viewing parameters changed that implies some
    //                  modifications of the buffers
    //                  (eg. type of normal, color/mono)
    //       * update() just to update the drawing

    // Call the base method to process others/classicals key
    Base::keyPressEvent(e);
  }

private:
  //! The window width in pixels.
  int m_width = CGAL_BASIC_VIEWER_INIT_SIZE_X;

  //! The window height in pixels.
  int m_height = CGAL_BASIC_VIEWER_INIT_SIZE_Y;

  //! The ratio between pixel and opengl units (in world coordinate system).
  double m_pixel_ratio = 1;

  //! The polygon with holes to draw.
  const Mpwh& m_mpwh;
};

// Specialization of draw function.
template<class T, class C>
void draw(const CGAL::Multipolygon_with_holes_2<T, C>& mpwh,
          const char* title = "Multipolygon_with_holes_2 Basic Viewer")
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite = true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (! cgal_test_suite) {
    using Mpwh = CGAL::Multipolygon_with_holes_2<T, C>;
    using Viewer = Mpwh_2_basic_viewer_qt<Mpwh>;
    CGAL::Qt::init_ogl_context(4,3);
    int argc = 1;
    const char* argv[2] = {"t2_viewer", nullptr};
    QApplication app(argc, const_cast<char**>(argv));
    Viewer mainwindow(app.activeWindow(), mpwh, title);
    mainwindow.add_elements();
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_MULTIPOLYGON_WITH_HOLES_2_H
