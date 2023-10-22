// Copyright(c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_DRAW_CT2_H
#define CGAL_DRAW_CT2_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_2/internal/In_domain.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/draw_triangulation_2.h>

namespace CGAL
{
// Specialization of draw function.
#define CGAL_T2_TYPE CGAL::Constrained_triangulation_2<Gt, Tds, Itag>

template <class Gt, class Tds, class Itag, class InDomainPmap, typename BufferType=float>
void add_in_graphics_scene(const CGAL_T2_TYPE& at2, InDomainPmap ipm,
                            CGAL::Graphics_scene& graphics_scene)
{
  using BASET2=CGAL::Triangulation_2<Gt, Tds>;

  Graphics_scene_options<BASET2, //CGAL_T2_TYPE,
                  typename CGAL_T2_TYPE::Vertex_handle,
                  typename CGAL_T2_TYPE::Finite_edges_iterator,
                  typename CGAL_T2_TYPE::Finite_faces_iterator>
    drawingFunctor;

  drawingFunctor.colored_edge =
    [](const BASET2& t2, typename CGAL_T2_TYPE::Finite_edges_iterator eh) -> bool
    { return static_cast<const CGAL_T2_TYPE&>(t2).is_constrained(*eh); };

  drawingFunctor.edge_color =
    [](const BASET2& t2, typename CGAL_T2_TYPE::Finite_edges_iterator eh) -> CGAL::IO::Color
  { return static_cast<const CGAL_T2_TYPE&>(t2).is_constrained(*eh)? CGAL::IO::green() : CGAL::IO::black(); };

  drawingFunctor.colored_face =
    [](const BASET2&, typename CGAL_T2_TYPE::Finite_faces_iterator) -> bool
    { return true; };

  drawingFunctor.face_color =
    [&ipm](const BASET2&, typename CGAL_T2_TYPE::Finite_faces_iterator fh) -> CGAL::IO::Color
  { return get(ipm, fh)? CGAL::IO::yellow() : CGAL::IO::white(); };

  add_in_graphics_scene(at2, graphics_scene, drawingFunctor);
}

template <class Gt, class Tds, class Itag, typename BufferType=float>
void add_in_graphics_scene(const CGAL_T2_TYPE& at2,
                            CGAL::Graphics_scene& graphics_scene)
{
  internal::In_domain<CGAL_T2_TYPE> in_domain;
  add_in_graphics_scene(at2, in_domain, graphics_scene);
}

template<class Gt, class Tds, class Itag, class InDomainPmap>
void draw(const CGAL_T2_TYPE& at2, InDomainPmap ipm,
          const char *title="Constrained Triangulation_2 Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_in_graphics_scene(at2, ipm, buffer);
  draw_graphics_scene(buffer, title);
}


template<class Gt, class Tds, class Itag>
void draw(const CGAL_T2_TYPE& at2,
          const char *title="Constrained Triangulation_2 Basic Viewer")
{
  internal::In_domain<CGAL_T2_TYPE> in_domain;
  draw(at2, in_domain, title);
}

#undef CGAL_T2_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_CT2_H
