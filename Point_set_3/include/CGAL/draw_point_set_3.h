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
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_DRAW_POINT_SET_3_H
#define CGAL_DRAW_POINT_SET_3_H

#include <CGAL/license/Point_set_3.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Point_set_3.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
\ingroup PkgDrawPointSet3D

opens a new window and draws `aps`, an instance of the `CGAL::Point_set_3` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam PS an instance of the `CGAL::Point_set_3` class.
\param aps the point set to draw.

*/
template<class PS>
void draw(const PS& aps);

} /* namespace CGAL */
#endif

namespace CGAL {

namespace draw_function_for_PointSet {

template <typename BufferType=float, class PointSet, class DrawingFunctor>
void compute_elements(const PointSet& pointset,
                      Graphic_buffer<BufferType>& graphic_buffer,
                      const DrawingFunctor& drawing_functor)
{
  if (!drawing_functor.are_vertices_enabled())
  { return; }

  for (typename PointSet::const_iterator it=pointset.begin();
       it!=pointset.end(); ++it)
  {
    if(drawing_functor.draw_vertex(pointset, it))
    {
      if (drawing_functor.colored_vertex(pointset, it))
      {
        graphic_buffer.add_point(pointset.point(*it),
                                 drawing_functor.vertex_color(pointset, it));
      }
      else
      { graphic_buffer.add_point(pointset.point(*it)); }
    }
  }
}

} // namespace draw_function_for_PointSet

template <class P, class V, typename BufferType=float, class DrawingFunctor>
void add_in_graphic_buffer(const Point_set_3<P, V>& apointset,
                           Graphic_buffer<BufferType>& graphic_buffer,
                           const DrawingFunctor& drawing_functor)
{
  draw_function_for_PointSet::compute_elements(apointset,
                                               graphic_buffer,
                                               drawing_functor);
}

template <class P, class V, typename BufferType=float>
void add_in_graphic_buffer(const Point_set_3<P, V>& apointset,
                           Graphic_buffer<BufferType>& graphic_buffer)
{
  CGAL::Drawing_functor<Point_set_3<P, V>,
                        typename Point_set_3<P, V>::const_iterator,
                        int, int> drawing_functor;
  add_in_graphic_buffer(apointset, graphic_buffer, drawing_functor);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function.
  template <class P, class V, class DrawingFunctor>
void draw(const Point_set_3<P, V>& apointset,
          const DrawingFunctor& drawing_functor,
          const char *title="Point_set_3 Basic Viewer")
{
  Graphic_buffer<float> buffer;
  add_in_graphic_buffer(apointset, buffer, drawing_functor);
  draw_buffer(buffer, title);
}

template <class P, class V>
void draw(const Point_set_3<P, V>& apointset,
          const char *title="Point_set_3 Basic Viewer")
{
  Graphic_buffer<float> buffer;
  add_in_graphic_buffer(apointset, buffer);
  draw_buffer(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

} // End namespace CGAL

#endif // CGAL_DRAW_POINT_SET_3_H
