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
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Random.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL {

namespace draw_function_for_FG {

template <typename BufferType=float, typename FG, typename DrawingFunctor>
void compute_elements(const FG &fg,
                      CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                      const DrawingFunctor &m_drawing_functor)
{
  using Point=typename boost::property_map_value<FG, CGAL::vertex_point_t>::type;
  using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Vector = typename Kernel::Vector_3;

  auto vnormals = get(CGAL::dynamic_vertex_property_t<Vector>(), fg);
  auto point_pmap = get(CGAL::vertex_point, fg);
  for (auto v : vertices(fg))
  {
    Vector n(NULL_VECTOR);
    int i = 0;
    for (auto h : halfedges_around_target(halfedge(v, fg), fg))
    {
      if (!is_border(h, fg))
      {
        Vector ni = CGAL::cross_product(
            Vector(get(point_pmap, source(h, fg)),
                    get(point_pmap, target(h, fg))),
            Vector(get(point_pmap, target(h, fg)),
                    get(point_pmap, target(next(h, fg), fg))));
        if (ni != NULL_VECTOR)
        {
          n += ni;
          ++i;
        }
      }
    }
    put(vnormals, v, n / i);
  }

  if (m_drawing_functor.are_faces_enabled())
  {
    for (auto fh : faces(fg))
    {
      if (fh != boost::graph_traits<FG>::null_face() && // face exists
          m_drawing_functor.colored_face && // std::function is not null
          m_drawing_functor.colored_face(fg, fh)) // and face is colored
      {
        graphic_buffer.face_begin(m_drawing_functor.face_color(fg, fh));
        auto hd = halfedge(fh, fg);
        const auto first_hd = hd;
        do
        {
          auto v = source(hd, fg);
          graphic_buffer.add_point_in_face(get(point_pmap, v), get(vnormals, v));
          hd = next(hd, fg);
        }
        while (hd != first_hd);
        graphic_buffer.face_end();
      }
    }
  }

  if(m_drawing_functor.are_edges_enabled())
  {
    for (auto e : edges(fg))
    {
      if(m_drawing_functor.colored_edge && // std::function is not null
         m_drawing_functor.colored_edge(fg, e)) // and edge is colored
      { 
        graphic_buffer.add_segment(get(point_pmap, source(halfedge(e, fg), fg)),
                                   get(point_pmap, target(halfedge(e, fg), fg)),
                                   m_drawing_functor.edge_color(fg, e));
      }
      else
      {
        graphic_buffer.add_segment(get(point_pmap, source(halfedge(e, fg), fg)),
                                   get(point_pmap, target(halfedge(e, fg), fg)));
      }
    }
  }

  if(m_drawing_functor.are_vertices_enabled())
  {
    for (auto v : vertices(fg))
    {
      if(m_drawing_functor.colored_vertex && // std::function is not null
         m_drawing_functor.colored_vertex(fg, v)) // and vertex is colored
      { 
        graphic_buffer.add_point(get(point_pmap, v),
                                 m_drawing_functor.vertex_color(fg, v));
      }
      else
      { 
        graphic_buffer.add_point(get(point_pmap, v));
      }
    }
  }
}

} // draw_function_for_FG

template <typename BufferType=float, class FG, class DrawingFunctor>
void add_in_graphic_buffer_for_fg(const FG &fg,
                                  CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                                  const DrawingFunctor &drawing_functor)
{
  draw_function_for_FG::compute_elements(fg, graphic_buffer, drawing_functor);
}

template <typename BufferType=float, class FG>
void add_in_graphic_buffer_for_fg(const FG &fg,
                                  CGAL::Graphic_buffer<BufferType> &graphic_buffer)
{
  Drawing_functor<FG,
                  typename boost::graph_traits<FG>::vertex_descriptor,
                  typename boost::graph_traits<FG>::edge_descriptor,
                  typename boost::graph_traits<FG>::face_descriptor>
      drawing_functor;

  drawing_functor.colored_face = [](const FG&,
             typename boost::graph_traits<FG>::face_descriptor) -> bool
  { return true; };


  drawing_functor.face_color =  [] (const FG&,
             typename boost::graph_traits<FG>::face_descriptor fh) -> CGAL::IO::Color
  {
    if (fh==boost::graph_traits<FG>::null_face())
    { return CGAL::IO::Color(100, 125, 200); }

    // TODO (?) use a seed given fh (cannot directly cast because FG is either a polyhedron or a surface mesh)
    return get_random_color(CGAL::get_default_random());
  };

  add_in_graphic_buffer_for_fg(fg, graphic_buffer, drawing_functor);
}

} // End namespace CGAL

#endif // CGAL_DRAW_SURFACE_MESH_H
