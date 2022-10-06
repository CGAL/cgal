// Copyright (c) 2018-2022 GeometryFactory (France)
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
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_DRAW_FACE_GRAPH_H
#define CGAL_DRAW_FACE_GRAPH_H

#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Random.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL {

namespace draw_function_for_SM {

template <typename BufferType = float, typename SM, typename DrawingFunctor>
void compute_elements(const SM &sm,
                      CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                      const DrawingFunctor &m_drawing_functor, bool anofaces = false) {

  using Point =
      typename boost::property_map_value<SM, CGAL::vertex_point_t>::type;
  using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Vector = typename Kernel::Vector_3;

  auto vnormals = get(CGAL::dynamic_vertex_property_t<Vector>(), sm);
  auto point_pmap = get(CGAL::vertex_point, sm);
  for (auto v : vertices(sm)) {
    Vector n(NULL_VECTOR);
    int i = 0;
    for (auto h : halfedges_around_target(halfedge(v, sm), sm)) {
      if (!is_border(h, sm)) {
        Vector ni = CGAL::cross_product(
            Vector(get(point_pmap, source(h, sm)),
                    get(point_pmap, target(h, sm))),
            Vector(get(point_pmap, target(h, sm)),
                    get(point_pmap, target(next(h, sm), sm))));
        if (ni != NULL_VECTOR) {
          n += ni;
          ++i;
        }
      }
    }
    put(vnormals, v, n / i);
  }

  if (!anofaces) {
    for (auto fh : faces(sm)) {
      if (fh != boost::graph_traits<SM>::null_face() && m_drawing_functor.colored_face(sm, fh)) {
        CGAL::IO::Color c = m_drawing_functor.face_color(sm, fh);
        graphic_buffer.face_begin(c);
        auto hd = halfedge(fh, sm);
        const auto first_hd = hd;
        do {
          auto v = source(hd, sm);
          graphic_buffer.add_point_in_face(get(point_pmap, v), get(vnormals, v));
          hd = next(hd, sm);
        } while (hd != first_hd);
        graphic_buffer.face_end();
      }
    }
  }

  for (auto e : edges(sm)) {
    graphic_buffer.add_segment(get(point_pmap, source(halfedge(e, sm), sm)),
                get(point_pmap, target(halfedge(e, sm), sm)));
  }

  for (auto v : vertices(sm)) {
    graphic_buffer.add_point(get(point_pmap, v));
  }
}

} // draw_function_for_SM

template <typename BufferType = float, class SM, class DrawingFunctor>
void add_in_graphic_buffer(const SM &sm, CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                             const DrawingFunctor &m_drawing_functor) {
  draw_function_for_SM::compute_elements(sm, graphic_buffer, m_drawing_functor);
}

template <typename BufferType = float, class SM>
void add_in_graphic_buffer(const SM &sm, CGAL::Graphic_buffer<BufferType> &graphic_buffer) {

  // Default functor; user can add his own functor.
  Drawing_functor<SM, typename boost::graph_traits<SM>::face_descriptor,
                  typename boost::graph_traits<SM>::face_descriptor,
                  typename boost::graph_traits<SM>::face_descriptor>
      drawing_functor;

  drawing_functor.colored_face = [](const SM &,
             typename boost::graph_traits<SM>::face_descriptor fh) -> bool
  { return true; };


  drawing_functor.face_color =  [] (const SM &,
             typename boost::graph_traits<SM>::face_descriptor fh) -> CGAL::IO::Color
  {
    if (fh ==
        boost::graph_traits<SM>::null_face()) // use to get the mono color
      return CGAL::IO::Color(100, 125, 200);     // R G B between 0-255

    return get_random_color(CGAL::get_default_random());
  };

  add_in_graphic_buffer(sm, graphic_buffer, drawing_functor);
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_SURFACE_MESH_H
