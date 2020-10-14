// Copyright (c) 2018-2020 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Laurent Rineau <laurent.rineau@cgal.org>

#ifndef CGAL_DRAW_FACE_GRAPH_H
#define CGAL_DRAW_FACE_GRAPH_H

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Random.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL
{

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorFaceGraph
{
  template<typename Graph>
  CGAL::Color operator()(const Graph&,
                         typename boost::graph_traits<Graph>::face_descriptor fh) const
  {
    if (fh==boost::graph_traits<Graph>::null_face()) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    return get_random_color(CGAL::get_default_random());
  }
};

class SimpleFaceGraphViewerQt : public Basic_viewer_qt
{
  using Base = Basic_viewer_qt;

public:
  SimpleFaceGraphViewerQt(QWidget* parent) :
    Base(parent, "", false, true, true, true, false),
    m_compute_elements_impl([]{})
  {
  }

  /// Construct the viewer.
  /// @param amesh the surface mesh to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  template <typename SM>
  SimpleFaceGraphViewerQt(QWidget* parent,
                          const SM& amesh,
                          const char* title="Basic Surface_mesh Viewer",
                          bool anofaces=false) :
    SimpleFaceGraphViewerQt(parent, amesh, title, anofaces, DefaultColorFunctorFaceGraph())
  {
  }

  template <typename SM, typename ColorFunctor>
  SimpleFaceGraphViewerQt(QWidget* parent,
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
  void set_face_graph(const SM& amesh,
                      bool anofaces,
                      ColorFunctor fcolor) {
    m_compute_elements_impl = compute_elements_functor(amesh, anofaces, fcolor);
  }

  template <typename SM, typename ColorFunctor>
  void set_face_graph(const SM& amesh,
                      bool anofaces=false) {
    set_mesh(amesh, anofaces, DefaultColorFunctorFaceGraph());
  }
protected:
  template <typename SM, typename ColorFunctor>
  std::function<void()>
  compute_elements_functor(const SM& sm,
                           bool anofaces,
                           ColorFunctor fcolor)
  {
    using Point =
        typename boost::property_map_value<SM, CGAL::vertex_point_t>::type;
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Vector = typename Kernel::Vector_3;

    auto vnormals = get(CGAL::dynamic_vertex_property_t<Vector>(), sm);
    auto point_pmap = get(CGAL::vertex_point, sm);
    for (auto v : vertices(sm))
    {
      Vector n(NULL_VECTOR);
      int i=0;
      for (auto h : halfedges_around_target(halfedge(v, sm), sm))
      {
        if (!is_border(h, sm))
        {
          Vector ni = CGAL::cross_product(
                        Vector(get(point_pmap, source(h, sm)), get(point_pmap, target(h, sm))),
                        Vector(get(point_pmap, target(h, sm)), get(point_pmap, target(next(h, sm), sm))));
          if (ni != NULL_VECTOR)
          {
            n+=ni;
            ++i;
          }
        }
      }
      put(vnormals, v, n/i);
    }

    // This function return a lambda expression, type-erased in a
    // `std::function<void()>` object.
    return [this, &sm, vnormals, anofaces, fcolor, point_pmap]()
    {
      this->clear();

      if (!anofaces)
      {
        for (auto fh: faces(sm))
        {
          if (fh!=boost::graph_traits<SM>::null_face())
          {
            CGAL::Color c=fcolor(sm, fh);
            face_begin(c);
            auto hd=halfedge(fh, sm);
            const auto first_hd = hd;
            do
              {
                auto v = source(hd, sm);
                add_point_in_face(get(point_pmap, v), get(vnormals, v));
                hd=next(hd, sm);
              }
            while(hd!=first_hd);
            face_end();
          }
        }
      }

      for (auto e: edges(sm))
      {
        add_segment(get(point_pmap, source(halfedge(e, sm), sm)),
                    get(point_pmap, target(halfedge(e, sm), sm)));
      }

      for (auto v: vertices(sm))
      {
        this->add_point(get(point_pmap, v));
      }
    };
  }

  void keyPressEvent(QKeyEvent *e) override
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

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_SURFACE_MESH_H
