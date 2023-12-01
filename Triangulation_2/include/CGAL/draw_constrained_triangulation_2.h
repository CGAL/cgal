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
#include <CGAL/draw_triangulation_2.h>

namespace CGAL
{

template<class CDT>
struct Graphics_scene_options_constrained_triangulation_2:
    public CGAL::Graphics_scene_options<typename CDT::Triangulation,
                                        typename CDT::Vertex_handle,
                                        typename CDT::Finite_edges_iterator,
                                        typename CDT::Finite_faces_iterator>
{
  using BASET2=typename CDT::Triangulation;

  Graphics_scene_options_constrained_triangulation_2(const CDT& cdt)
  {
    this->colored_edge =
      [&cdt](const BASET2&, typename CDT::Finite_edges_iterator eh) -> bool
      { return cdt.is_constrained(*eh); };

    this->edge_color =
      [&cdt](const BASET2&, typename CDT::Finite_edges_iterator eh) -> CGAL::IO::Color
      { return cdt.is_constrained(*eh)? CGAL::IO::green() : CGAL::IO::black(); };
  };

  template<class InDomainPmap>
  Graphics_scene_options_constrained_triangulation_2(const CDT& cdt, InDomainPmap ipm)
  {
    this->colored_edge =
      [&cdt](const BASET2&, typename CDT::Finite_edges_iterator eh) -> bool
      { return cdt.is_constrained(*eh); };

    this->edge_color =
      [&cdt](const BASET2&, typename CDT::Finite_edges_iterator eh) -> CGAL::IO::Color
      { return cdt.is_constrained(*eh)? CGAL::IO::green() : CGAL::IO::black(); };

    this->colored_face =
      [](const BASET2&, typename CDT::Finite_faces_iterator) -> bool
      { return true; };

    this->face_color =
      [ipm](const BASET2&, typename CDT::Finite_faces_iterator fh) -> CGAL::IO::Color
      { return get(ipm, fh)? CGAL::IO::blue() : CGAL::IO::white(); };
  };
};

// Specialization of draw function.
#define CGAL_T2_TYPE CGAL::Constrained_triangulation_2<Gt, Tds, Itag>

template <class Gt, class Tds, class Itag, class InDomainPmap>
void add_to_graphics_scene(const CGAL_T2_TYPE& at2, InDomainPmap ipm,
                           CGAL::Graphics_scene& graphics_scene)
{
  Graphics_scene_options_constrained_triangulation_2<CGAL_T2_TYPE> gso(at2, ipm);
  draw_function_for_t2::compute_elements(at2, graphics_scene, gso);
}

template <class Gt, class Tds, class Itag>
void add_to_graphics_scene(const CGAL_T2_TYPE& at2,
                           CGAL::Graphics_scene& graphics_scene)
{
  Graphics_scene_options_constrained_triangulation_2<CGAL_T2_TYPE> gso(at2);
  draw_function_for_t2::compute_elements(at2, graphics_scene, gso);
}

#ifdef CGAL_USE_BASIC_VIEWER

template<class Gt, class Tds, class Itag, class InDomainPmap>
void draw(const CGAL_T2_TYPE& at2, InDomainPmap ipm,
          const char *title="Constrained Triangulation_2 Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(at2, ipm, buffer);
  draw_graphics_scene(buffer, title);
}

template<class Gt, class Tds, class Itag>
void draw(const CGAL_T2_TYPE& at2,
          const char *title="Constrained Triangulation_2 Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(at2, buffer);
  draw_graphics_scene(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_T2_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_CT2_H
