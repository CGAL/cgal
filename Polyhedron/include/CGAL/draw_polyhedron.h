// Copyright (c) 2018-2020  ETH Zurich (Switzerland).
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

#ifndef CGAL_DRAW_POLYHEDRON_H
#define CGAL_DRAW_POLYHEDRON_H

#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/license/Polyhedron.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_face_graph.h>
#include <CGAL/Random.h>

namespace CGAL
{

// Specialization of draw function.
#define CGAL_POLY_TYPE CGAL::Polyhedron_3 \
  <PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>

template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc, typename BufferType = float>
void draw(const CGAL_POLY_TYPE& apoly,
          const char* title="Polyhedron Basic Viewer",
          bool nofill=false)
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(apoly, buffer);
  draw_buffer(buffer);
}

template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc, typename BufferType = float, class DrawingFunctor>
void draw(const CGAL_POLY_TYPE& apoly,
          const DrawingFunctor &drawing_functor,
          const char* title="Polyhedron Basic Viewer",
          bool nofill=false)
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(apoly, buffer, drawing_functor);
  draw_buffer(buffer);
}

#undef CGAL_POLY_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYHEDRON_H
