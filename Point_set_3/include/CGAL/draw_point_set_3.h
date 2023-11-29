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
#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Point_set_3.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
\ingroup PkgDrawPointSet3D

opens a new window and draws `aps`, an instance of the `CGAL::Point_set_3` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam PS an instance of the `CGAL::Point_set_3` class.
\param aps the point set to draw.

*/
template<class PS>
void draw(const PS& aps);

} /* namespace CGAL */
#endif

namespace CGAL {

namespace draw_function_for_PointSet {

template <class PointSet, class GSOptions>
void compute_elements(const PointSet& pointset,
                      Graphics_scene& graphics_scene,
                      const GSOptions& gs_options)
{
  if (!gs_options.are_vertices_enabled())
  { return; }

  for (typename PointSet::const_iterator it=pointset.begin();
       it!=pointset.end(); ++it)
  {
    if(gs_options.draw_vertex(pointset, it))
    {
      if (gs_options.colored_vertex(pointset, it))
      {
        graphics_scene.add_point(pointset.point(*it),
                                 gs_options.vertex_color(pointset, it));
      }
      else
      { graphics_scene.add_point(pointset.point(*it)); }
    }
  }
}

} // namespace draw_function_for_PointSet

template <class P, class V, class GSOptions>
void add_to_graphics_scene(const Point_set_3<P, V>& apointset,
                            Graphics_scene& graphics_scene,
                            const GSOptions& gs_options)
{
  draw_function_for_PointSet::compute_elements(apointset,
                                               graphics_scene,
                                               gs_options);
}

template <class P, class V>
void add_to_graphics_scene(const Point_set_3<P, V>& apointset,
                            Graphics_scene& graphics_scene)
{
  CGAL::Graphics_scene_options<Point_set_3<P, V>,
                               typename Point_set_3<P, V>::const_iterator,
                               int, int> gs_options;
  add_to_graphics_scene(apointset, graphics_scene, gs_options);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function.
  template <class P, class V, class GSOptions>
void draw(const Point_set_3<P, V>& apointset,
          const GSOptions& gs_options,
          const char *title="Point_set_3 Basic Viewer")
{
  Graphics_scene buffer;
  add_to_graphics_scene(apointset, buffer, gs_options);
  draw_graphics_scene(buffer, title);
}

template <class P, class V>
void draw(const Point_set_3<P, V>& apointset,
          const char *title="Point_set_3 Basic Viewer")
{
  Graphics_scene buffer;
  add_to_graphics_scene(apointset, buffer);
  draw_graphics_scene(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

} // End namespace CGAL

#endif // CGAL_DRAW_POINT_SET_3_H
