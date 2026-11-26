// Copyright (c) 2012 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel          <efifogel@gmail.com>
//            Shepard Liu        <shepard0liu@gmail.com>

#ifndef CGAL_DRAW_ARRANGEMENT_2_H
#define CGAL_DRAW_ARRANGEMENT_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <type_traits>

#ifdef CGAL_USE_BASIC_VIEWER
#include <QApplication>
#endif

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Basic_viewer.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/assertions.h>
#include <CGAL/config.h>
#include <CGAL/Draw_aos/Arr_graphics_scene_builder.h>
#include <CGAL/Draw_aos/type_utils.h>

#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/Draw_aos/Arr_viewer.h>
#endif

namespace CGAL {
namespace draw_aos {

// Detect whether T is or derives from Arrangement_on_surface_2<*, *>
template <typename GeometryTraits_2, typename TopologyTraits>
std::true_type is_aos_test(const Arrangement_on_surface_2<GeometryTraits_2, TopologyTraits>*);

std::false_type is_aos_test(...);

template <typename T>
inline constexpr bool is_or_derived_from_aos_v = decltype(is_aos_test(static_cast<const T*>(nullptr)))::value;

template <typename T>
using enable_if_aos_t = std::enable_if_t<is_or_derived_from_aos_v<T>, int>;

//!
inline void draw_unimplemented() {
  CGAL_error_msg("Geometry traits type of arrangement is required to support approximation of Point_2 and "
                 "X_monotone_curve_2. Traits on curved surfaces needs additional support for parameterization.");
}

//! determines whether a parameter space bounding box is degenerate.
inline bool is_degenerate(const Bbox_2& bbox) { return (bbox.x_span() <= 0) || (bbox.y_span() <= 0); }

} // namespace draw_aos

/*! draws an arrangement on surface.
 *
 * The viewer approximates the arrangement within the visible part of the
 * parameter space and re-approximates it whenever the camera changes.
 *
 * \tparam Arrangement an instance of `Arrangement_on_surface_2` or a type derived from such an instance.
 * \tparam GSOptions a model of the `GraphicsSceneOptions` concept.
 * \param arr the arrangement to be drawn
 * \param gso graphics scene options
 * \param title title of the viewer window
 * \param initial_bbox parameter space bounding box to be shown initially. If empty, the approximate bounding box of the
 * arrangement is used. For arrangements induced by unbounded curves, the default initial bounding box is computed from
 * the bounded parts of the curves.
 */
template <typename Arrangement, typename GSOptions, draw_aos::enable_if_aos_t<Arrangement> = 0>
void draw(const Arrangement& arr,
          const GSOptions& gso,
          const char* title = "2D Arrangement on Surface Viewer",
          Bbox_2 initial_bbox = Bbox_2(0, 0, 0, 0)) {
  CGAL_USE(arr);
  CGAL_USE(gso);
  CGAL_USE(title);
  CGAL_USE(initial_bbox);

#ifdef CGAL_USE_BASIC_VIEWER
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  if constexpr (! draw_aos::has_approximate_traits_v<Geom_traits>) {
    draw_aos::draw_unimplemented();
  }
  else {
#if defined(CGAL_TEST_SUITE)
    bool cgal_test_suite = true;
#else
    bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif
    if (cgal_test_suite) return;

    Qt::init_ogl_context(4, 3);
    int argc = 1;
    const char* argv[2] = {title, nullptr};
    QApplication app(argc, const_cast<char**>(argv));
    draw_aos::Arr_viewer<Arrangement, GSOptions> viewer(app.activeWindow(), arr, gso, title, initial_bbox);
    viewer.show();
    app.exec();
  }
#endif
}

/*! draws an arrangement on surface with default graphics scene options. Faces are colored randomly.
 *
 * \tparam Arrangement an instance of `Arrangement_on_surface_2` or a type derived from such an instance.
 * \param arr the arrangement to be drawn
 * \param title title of the viewer window
 * \param initial_bbox parameter space bounding box to be shown initially. If empty, the approximate bounding box of the
 * arrangement is used. For arrangements induced by unbounded curves, the default initial bounding box is computed from
 * the bounded parts of the curves.
 */
template <typename Arrangement, draw_aos::enable_if_aos_t<Arrangement> = 0>
void draw(const Arrangement& arr,
          const char* title = "2D Arrangement on Surface Viewer",
          Bbox_2 initial_bbox = Bbox_2(0, 0, 0, 0))
{ draw(arr, draw_aos::default_gs_options<Arrangement>(), title, initial_bbox); }

/*! adds the vertices, edges, and faces of an arrangement on surface to a
 * graphics scene. The scene is not cleared before it is filled. The faces are
 * triangulated.
 *
 * As opposed to `draw()`, which re-approximates the arrangement whenever the
 * camera changes, the arrangement is approximated once. The approximation
 * error is set to a small fraction of the bounding box.
 *
 * \tparam Arrangement an instance of `Arrangement_on_surface_2` or a type derived from such an instance.
 * \tparam GSOptions a model of the `GraphicsSceneOptions` concept.
 * \param arr the arrangement to add to the graphics scene.
 * \param gs the graphics scene to fill.
 * \param gso the graphics scene options.
 * \param bbox the parameter space bounding box to which the arrangement is clipped. If empty, the approximate bounding
 * box of the arrangement is used. For arrangements induced by unbounded curves, the default bounding box is computed
 * from the bounded parts of the curves. Ignored for arrangements on curved surfaces, which are added entirely.
 */
template <typename Arrangement, typename GSOptions, draw_aos::enable_if_aos_t<Arrangement> = 0>
void add_to_graphics_scene(const Arrangement& arr,
                           Graphics_scene& gs,
                           const GSOptions& gso,
                           const Bbox_2& bbox = Bbox_2(0, 0, 0, 0)) {
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  if constexpr (! draw_aos::has_approximate_traits_v<Geom_traits>) {
    CGAL_USE(arr);
    CGAL_USE(gs);
    CGAL_USE(gso);
    CGAL_USE(bbox);
    draw_aos::draw_unimplemented();
  }
  else {
    draw_aos::Arr_graphics_scene_builder<Arrangement, GSOptions> builder(arr, gso);
    constexpr bool is_curved = draw_aos::is_or_derived_from_curved_surf_traits_v<Geom_traits>;
    Bbox_2 bb = (is_curved || draw_aos::is_degenerate(bbox)) ? builder.bbox() : bbox;
    builder(gs, bb, builder.approximation_error(bb));
  }
}

/*! adds the vertices, edges, and faces of an arrangement on surface to a
 * graphics scene with default graphics scene options. Faces are colored
 * randomly.
 *
 * \tparam Arrangement an instance of `Arrangement_on_surface_2` or a type derived from such an instance.
 * \param arr the arrangement to add to the graphics scene.
 * \param gs the graphics scene to fill.
 * \param bbox the parameter space bounding box to which the arrangement is clipped; see above.
 */
template <typename Arrangement, draw_aos::enable_if_aos_t<Arrangement> = 0>
void add_to_graphics_scene(const Arrangement& arr, Graphics_scene& gs, const Bbox_2& bbox = Bbox_2(0, 0, 0, 0))
{ add_to_graphics_scene(arr, gs, draw_aos::default_gs_options<Arrangement>(), bbox); }

} // namespace CGAL

#endif
