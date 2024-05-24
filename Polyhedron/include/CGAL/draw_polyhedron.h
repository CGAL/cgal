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

#include <CGAL/license/Polyhedron.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_face_graph.h>
#include <CGAL/Qt/Basic_viewer.h>

namespace CGAL
{

#define CGAL_POLY_TYPE CGAL::Polyhedron_3 \
  <PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>

// Specialization of add_to_graphics_scene function.
template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc,
         class GSOptions>
void add_to_graphics_scene(const CGAL_POLY_TYPE& apoly,
                           CGAL::Graphics_scene &graphics_scene,
                           const GSOptions &gs_options)
{ add_to_graphics_scene_for_fg(apoly, graphics_scene, gs_options); }

template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc>
void add_to_graphics_scene(const CGAL_POLY_TYPE& apoly,
                           CGAL::Graphics_scene &graphics_scene)
{ add_to_graphics_scene_for_fg(apoly, graphics_scene); }

// Specialization of draw function: require Qt and the CGAL basic viewer.
#ifdef CGAL_USE_BASIC_VIEWER

template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc>
void draw(const CGAL_POLY_TYPE& apoly,
          const char* title="Polyhedron Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene_for_fg(apoly, buffer);
  draw_graphics_scene(buffer, title);
}

template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc,
         class GSOptions>
void draw(const CGAL_POLY_TYPE& apoly,
          const GSOptions &gs_options,
          const char* title="Polyhedron Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene_for_fg(apoly, buffer, gs_options);
  draw_graphics_scene(buffer, title);
}
#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_POLY_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_POLYHEDRON_H
