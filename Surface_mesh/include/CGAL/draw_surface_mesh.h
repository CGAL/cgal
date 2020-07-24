// Copyright (c) 2018 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_SURFACE_MESH_H
#define CGAL_DRAW_SURFACE_MESH_H

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Dynamic_property_map.h>

#ifdef DOXYGEN_RUNNING

/*!
\ingroup PkgDrawSurfaceMesh

Open a new window and draw `asm`, an instance of the `CGAL::Surface_mesh` class. The function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam SM an instance of the `CGAL::Surface_mesh` class.
\param asm the surface mesh to draw.

*/
template<class SM>
void draw(const SM& asm);

#else // DOXYGEN_RUNNING

#include <CGAL/license/Surface_mesh.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Surface_mesh.h>
#include <CGAL/Random.h>

namespace CGAL
{

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorSM
{
  template<typename SM>
  CGAL::Color operator()(const SM&,
                         typename SM::Face_index fh) const
  {
    if (fh==boost::graph_traits<SM>::null_face()) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    CGAL::Random random((unsigned int)fh);
    return get_random_color(random);
  }
};

class SimpleSurfaceMeshViewerQt : public Basic_viewer_qt
{
  using Base = Basic_viewer_qt;

public:
  /// Construct the viewer.
  /// @param amesh the surface mesh to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  template <typename SM>
  SimpleSurfaceMeshViewerQt(QWidget* parent,
                            const SM& amesh,
                            const char* title="Basic Surface_mesh Viewer",
                            bool anofaces=false)
    : SimpleSurfaceMeshViewerQt(parent, amesh, title, anofaces, DefaultColorFunctorSM())
  {
  }

  template <typename SM, typename ColorFunctor>
  SimpleSurfaceMeshViewerQt(QWidget* parent,
                            const SM& amesh,
                            const char* title,
                            bool anofaces,
                            ColorFunctor fcolor) :
    // First draw: no vertex; edges, faces; mono-color; inverse normal
    Base(parent, title, false, true, true, true, false),
    m_compute_elements_impl(compute_elements_functor(amesh, anofaces, fcolor))
  {
  }

  void init() override {
    compute_elements();
    Base::init();
  }

  void compute_elements() {
    m_compute_elements_impl();
  }

  template <typename SM, typename ColorFunctor>
  void set_mesh(const SM& amesh,
                bool anofaces=false,
                ColorFunctor fcolor=DefaultColorFunctorSM()) {
    m_compute_elements_impl = compute_elements_functor(amesh, anofaces, fcolor);
    redraw();
  }

protected:
  template <typename SM, typename ColorFunctor>
  std::function<void()>
  compute_elements_functor(const SM& sm,
                           bool anofaces,
                           ColorFunctor fcolor)
  {
    using Point = typename SM::Point;
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Vector = typename Kernel::Vector_3;

    auto vnormals = get(CGAL::dynamic_vertex_property_t<Vector>(), sm);
    {
      // temporary face property map needed by `compute_normals`
      auto fpm = get(CGAL::dynamic_face_property_t<Vector>(), sm);

      CGAL::Polygon_mesh_processing::compute_normals(sm, vnormals, fpm);
    }

    // This function return a lambda expression, type-erased in a
    // `std::function<void()>` object.
    return [this, &sm, vnormals, anofaces, fcolor]()
    {
      this->clear();

      if (!anofaces)
      {
        for (typename SM::Face_range::iterator f=sm.faces().begin();
             f!=sm.faces().end(); ++f)
        {
          if (*f!=boost::graph_traits<SM>::null_face())
            { this->compute_face(sm, vnormals, *f, fcolor); }
        }
      }

      for (typename SM::Edge_range::iterator e=sm.edges().begin();
           e!=sm.edges().end(); ++e)
      { this->compute_edge(sm, *e); }

      for (typename SM::Vertex_range::iterator v=sm.vertices().begin();
           v!=sm.vertices().end(); ++v)
      { this->compute_vertex(sm, *v); }
    };
  }

  template <typename SM, typename VNormals, typename face_descriptor, typename ColorFunctor>
  void compute_face(const SM& sm, VNormals vnormals,
                    face_descriptor fh, const ColorFunctor& fcolor)
  {
    CGAL::Color c=fcolor(sm, fh);
    face_begin(c);
    auto hd=sm.halfedge(fh);
    do
    {
      auto v = sm.source(hd);
      add_point_in_face(sm.point(v), get(vnormals, v));
      hd=sm.next(hd);
    }
    while(hd!=sm.halfedge(fh));
    face_end();
  }

  template <typename SM, typename edge_descriptor>
  void compute_edge(const SM& sm, edge_descriptor e)
  {
    add_segment(sm.point(sm.source(sm.halfedge(e))),
                sm.point(sm.target(sm.halfedge(e))));
  }

  template <typename SM, typename vertex_descriptor>
  void compute_vertex(const SM& sm, vertex_descriptor vh)
  { add_point(sm.point(vh)); }


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
  std::function<void()> m_compute_elements_impl;
};

// Specialization of draw function.
template<class K>
void draw(const Surface_mesh<K>& amesh,
          const char* title="Surface_mesh Basic Viewer",
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
    const char* argv[2]={"surface_mesh_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    SimpleSurfaceMeshViewerQt mainwindow(app.activeWindow(), amesh, title,
                                         nofill);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // DOXYGEN_RUNNING

#endif // CGAL_DRAW_SURFACE_MESH_H
