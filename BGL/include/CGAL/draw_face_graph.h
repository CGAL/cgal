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

namespace CGAL {

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorFaceGraph
{
  template<typename Graph>
  CGAL::IO::Color operator()(const Graph& /*g*/,
                             typename boost::graph_traits<Graph>::face_descriptor /*f*/) const
  {
    return get_random_color(CGAL::get_default_random());
  }

  // edges and vertices are black by default
  template<typename Graph>
  CGAL::IO::Color operator()(const Graph& /*g*/,
                             typename boost::graph_traits<Graph>::edge_descriptor /*e*/) const
  {
    return IO::black();
  }

  template<typename Graph>
  CGAL::IO::Color operator()(const Graph& /*g*/,
                             typename boost::graph_traits<Graph>::vertex_descriptor /*v*/) const
  {
    return IO::black();
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
  /// @param g the face graph to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        useful for very big objects where this time could be long)
  template <typename Graph>
  SimpleFaceGraphViewerQt(QWidget* parent,
                          const Graph& g,
                          const char* title="Basic Face Graph Viewer",
                          bool anofaces=false) :
    SimpleFaceGraphViewerQt(parent, g, title, anofaces, DefaultColorFunctorFaceGraph())
  {
  }

  template <typename Graph, typename ColorFunctor>
  SimpleFaceGraphViewerQt(QWidget* parent,
                          const Graph& g,
                          const char* title,
                          bool anofaces,
                          ColorFunctor fcolor) :
    // First draw: no vertex; edges, faces; mono-color; inverse normal
    Base(parent, title, false, true, true, true, false),
    m_compute_elements_impl(compute_elements_functor(g, anofaces, fcolor))
  {
  }

  void init() override {
    compute_elements();
    Base::init();
  }

  void compute_elements() {
    m_compute_elements_impl();
  }

  template <typename Graph, typename ColorFunctor>
  void set_face_graph(const Graph& g,
                      bool anofaces,
                      ColorFunctor fcolor) {
    m_compute_elements_impl = compute_elements_functor(g, anofaces, fcolor);
  }

  template <typename Graph, typename ColorFunctor>
  void set_face_graph(const Graph& g,
                      bool anofaces=false) {
    set_mesh(g, anofaces, DefaultColorFunctorFaceGraph());
  }
protected:
  template <typename Graph, typename ColorFunctor>
  std::function<void()>
  compute_elements_functor(const Graph& g,
                           bool anofaces,
                           ColorFunctor fcolor)
  {
    using Point = typename boost::property_map_value<Graph, CGAL::vertex_point_t>::type;
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Vector = typename Kernel::Vector_3;

    auto vnormals = get(CGAL::dynamic_vertex_property_t<Vector>(), g);
    auto point_pmap = get(CGAL::vertex_point, g);
    for (auto v : vertices(g))
    {
      Vector n(NULL_VECTOR);
      int i=0;
      for (auto h : halfedges_around_target(halfedge(v, g), g))
      {
        if (!is_border(h, g))
        {
          Vector ni = CGAL::cross_product(
                        Vector(get(point_pmap, source(h, g)), get(point_pmap, target(h, g))),
                        Vector(get(point_pmap, target(h, g)), get(point_pmap, target(next(h, g), g))));
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
    return [this, &g, vnormals, anofaces, fcolor, point_pmap]()
    {
      this->clear();

      if (!anofaces)
      {
        for (auto fh: faces(g))
        {
          const CGAL::IO::Color& c = fcolor(g, fh);
          face_begin(c);
          auto hd=halfedge(fh, g);
          const auto first_hd = hd;
          do
            {
              auto v = source(hd, g);
              add_point_in_face(get(point_pmap, v), get(vnormals, v));
              hd=next(hd, g);
            }
          while(hd!=first_hd);
          face_end();
        }
      }

      for (auto e: edges(g))
      {
        const CGAL::IO::Color& c = fcolor(g, e);
        add_segment(get(point_pmap, source(halfedge(e, g), g)),
                    get(point_pmap, target(halfedge(e, g), g)),
                    c);
      }

      for (auto v: vertices(g))
      {
        const CGAL::IO::Color& c = fcolor(g, v);
        this->add_point(get(point_pmap, v), c);
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
