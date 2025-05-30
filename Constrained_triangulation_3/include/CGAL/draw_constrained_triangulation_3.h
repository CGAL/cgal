// Copyright (c) 2025 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_DRAW_CONSTRAINED_T3_H
#define CGAL_DRAW_CONSTRAINED_T3_H

#include <CGAL/license/Constrained_triangulation_3.h>
#include <CGAL/draw_triangulation_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_3_fwd.h>
#include <CGAL/type_traits.h>
#include <algorithm>

namespace CGAL {
/*!
\ingroup PkgDrawCDT_3

opens a new window and draws the constrained triangulation.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.
*/
template <typename Traits, typename Tr>
void draw(const Conforming_constrained_Delaunay_triangulation_3<Traits, Tr>& ccdt,
          const char *title="Constrained 3D Triangulation Basic Viewer")
{
  using Tr_ = CGAL::cpp20::remove_cvref_t<decltype(ccdt.triangulation())>;
  using Vertex_handle = typename Tr_::Vertex_handle;
  using Cell_handle = typename Tr_::Cell_handle;
  using Edge_descriptor = typename Tr_::Finite_edges_iterator;
  using Facet_descriptor = typename Tr_::Finite_facets_iterator;

  using Face_index = CGAL::cpp20::remove_cvref_t<
      decltype(std::declval<Cell_handle>()->ccdt_3_data().face_constraint_index(0))>;

  Face_index nb_colors = 0;
  std::for_each(
      ccdt.constrained_facets_begin(), ccdt.constrained_facets_end(),
      [&](const auto& f) {
        auto [c, index] = f;
        nb_colors = (std::max)(nb_colors, c->ccdt_3_data().face_constraint_index(index) + 1);
      });
  std::vector<CGAL::IO::Color> colors(nb_colors);
  std::generate(colors.begin(), colors.end(), []() {
    return CGAL::get_random_color(CGAL::get_default_random());
  });
  CGAL::Graphics_scene_options<Tr_, Vertex_handle, Edge_descriptor, Facet_descriptor> options;
  options.draw_face = [](const Tr_&, Facet_descriptor f) {
    auto [c, index] = *f;
    return c->ccdt_3_data().is_facet_constrained(index);
  };
  options.colored_face = [](const Tr_&, Facet_descriptor) {
    return true;
  };
  options.face_color = [&](const Tr_&, Facet_descriptor f) {
    auto [c, index] = *f;
    return colors[c->ccdt_3_data().face_constraint_index(index)];
  };
  draw(ccdt.triangulation(), options, title);
}

/*!
\ingroup PkgDrawCDT_3

A shortcut to \link PkgDrawTriangulation3 `CGAL::draw(ccdt.triangulation(), gs_options, title)` \endlink.
*/
template <typename Traits, typename Tr, typename GSOptions>
void draw(const Conforming_constrained_Delaunay_triangulation_3<Traits, Tr>& ccdt,
          const GSOptions& gs_options,
          const char *title="Constrained 3D Triangulation Basic Viewer")
{
  draw(ccdt.triangulation(), gs_options, title);
}

} // End namespace CGAL

#endif // CGAL_DRAW_CONSTRAINED_T3_H
