// Copyright (c) 2018-2020 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_DRAW_SURFACE_MESH_H
#define CGAL_DRAW_SURFACE_MESH_H

#ifdef DOXYGEN_RUNNING

/*!
\ingroup PkgDrawSurfaceMesh

Open a new window and draw `asm`, an instance of the `CGAL::Surface_mesh` class. The function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam SM an instance of the `CGAL::Surface_mesh` class.
\param asm the surface mesh to draw.

*/
template<class SM>
void draw(const SM& asm);

#else // DOXYGEN_RUNNING

#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/license/Surface_mesh.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_face_graph.h>

namespace CGAL
{

// Specialization of draw function.
template<class K, typename BufferType = float>
void draw(const Surface_mesh<K>& amesh,
          const char* title="Surface_mesh Basic Viewer",
          bool nofill=false)
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(amesh, buffer);
  draw_buffer(buffer);
}

template<class K, typename BufferType = float, class DrawingFunctor>
void draw(const Surface_mesh<K>& amesh,
          const DrawingFunctor &drawing_functor,
          const char* title="Surface_mesh Basic Viewer",
          bool nofill=false)
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(amesh, buffer, drawing_functor);
  draw_buffer(buffer);
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // DOXYGEN_RUNNING

#endif // CGAL_DRAW_SURFACE_MESH_H
