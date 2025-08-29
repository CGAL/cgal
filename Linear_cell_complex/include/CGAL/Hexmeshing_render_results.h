// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>, Théo Bénard <benard320@gmail.com>
//
#ifndef HEXMESHING_RENDER_RESULTS_H
#define HEXMESHING_RENDER_RESULTS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/IO/Color.h>


namespace CGAL {
  template <typename LCC>
  using LCCSceneOptions = CGAL::Graphics_scene_options<LCC,
                          typename LCC::Dart_const_handle,
                          typename LCC::Dart_const_handle,
                          typename LCC::Dart_const_handle,
                          typename LCC::Dart_const_handle>;

  inline IO::Color viridis255(double t) {
    t = std::clamp(t, 0.0, 1.0);
    static const double c[5][3] = {
      {68, 1, 84}, // 0.00
      {59, 82, 139}, // 0.25
      {33, 145, 140}, // 0.50
      {94, 201, 98}, // 0.75
      {253, 231, 37}  // 1.00
    };
    const double pos = t * 4.0;
    const int i      = int(pos);
    const double s   = pos - i;
    const auto lerp  = [&](int k){ return (1-s)*c[i][k] + s*c[i+1][k]; };
    return { uint8_t(lerp(0)), uint8_t(lerp(1)), uint8_t(lerp(2)) };
  }

  /**
   * @brief Renders the result of the two-refinement algorithm as a graphics scene
   * 
   * This function creates a visual representation of the hexahedral mesh generated
   * by the two-refinement algorithm. 
   * 
   * The function performs the following operations:
   * 
   * 1. **Scene Options Setup**: Creates `LCCSceneOptions` with custom rendering behavior:
   *    - Enables colored volume rendering for all volumes
   *    - Assigns random colors to volumes using `rand_color_from_dart`
   *    - Configures volume visibility based on intersection with AABB tree when trimming is enabled
   * 
   * 2. **Graphics Scene Creation**: Creates a `CGAL::Graphics_scene` buffer and adds
   *    the Linear Cell Complex to it using the configured scene options
   * 
   * 3. **Scene Display**: Calls `CGAL::draw_graphics_scene` to render the mesh
   * 
   * @param hdata The hex-meshing data which contains refined linear cell complex
   * @param title Title for the graphics window (default: "TwoRefinement Result")
   */
  void render_two_refinement_result(const Hexmeshing_for_linear_cell_complex& hdata, const char* title = "TwoRefinement Result"){
    const internal::Hexmeshing::LCC& lcc = hdata.lcc;

    LCCSceneOptions<internal::Hexmeshing::LCC> gso;

    gso.colored_volume = [&](const internal::Hexmeshing::LCC& lcc, internal::Hexmeshing::LCC::Dart_const_handle dart){ return true; };
    gso.volume_color = [&](const internal::Hexmeshing::LCC& lcc, internal::Hexmeshing::LCC::Dart_const_handle dart){
      double f = lcc.attribute<3>(dart)->info().fraction;
      return viridis255(f);
    };
    gso.draw_volume = [&](const internal::Hexmeshing::LCC& lcc, internal::Hexmeshing::LCC::Dart_const_handle dart){
      return true;
    };

    Graphics_scene buffer;
    add_to_graphics_scene(lcc, buffer, gso);
    draw_graphics_scene(buffer);
  }

  void render_two_refinement_result_with_mark(const Hexmeshing_for_linear_cell_complex& hdata, internal::Hexmeshing::size_type mark, const char* title = "TwoRefinement Result"){
    const internal::Hexmeshing::LCC& lcc = hdata.lcc;
    
    LCCSceneOptions<internal::Hexmeshing::LCC> gso;

    gso.colored_volume = [&](const internal::Hexmeshing::LCC& lcc, internal::Hexmeshing::LCC::Dart_const_handle dart){ return true; };
    gso.volume_color = [&](const internal::Hexmeshing::LCC& lcc, internal::Hexmeshing::LCC::Dart_const_handle dart){
      return viridis255(lcc.attribute<3>(dart)->info().fraction);
    };
    gso.draw_volume = [&](const internal::Hexmeshing::LCC& lcc, internal::Hexmeshing::LCC::Dart_const_handle dart){
      return lcc.is_marked(dart, mark);
    };

    Graphics_scene buffer;
    add_to_graphics_scene(lcc, buffer, gso);
    draw_graphics_scene(buffer);
  }
}



#endif