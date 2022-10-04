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

#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/license/Point_set_3.h>

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

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Point_set_3.h>
#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Random.h>

namespace CGAL {

namespace draw_function_for_PointSet {

template <typename BufferType = float, class PointSet>
void compute_vertex(const typename PointSet::Point_map::value_type &p,
                    Graphic_buffer<BufferType> &graphic_buffer) {
  graphic_buffer.add_point(p);

  // We can use add_point(p, c) with c a CGAL::IO::Color to add a colored point
  // E.g: graphic_buffer.add_point(p, CGAL::IO::Color(100, 125, 200));
}

template <typename BufferType = float, class PointSet>
void compute_elements(const PointSet *pointset,
                      Graphic_buffer<BufferType> &graphic_buffer) {
  for (typename PointSet::const_iterator it = pointset->begin();
       it != pointset->end(); ++it) {
    compute_vertex<float, PointSet>(pointset->point(*it), graphic_buffer);
  }
}

} // namespace draw_function_for_PointSet

template <typename BufferType = float, class PointSet>
void add_in_graphic_buffer(Graphic_buffer<BufferType> &graphic_buffer,
                                     const PointSet *aPointSet = nullptr) {
  if (aPointSet != nullptr) {
    draw_function_for_PointSet::compute_elements(aPointSet, graphic_buffer);
  }
}

// Specialization of draw function.
template <class P, class V>
void draw(const Point_set_3<P, V> &apointset,
          const char *title = "Point_set_3 Basic Viewer") {
  Graphic_buffer<float> buffer;
  add_in_graphic_buffer(buffer, &apointset);
  draw_buffer(buffer);
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POINT_SET_3_H
