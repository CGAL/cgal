// Copyright (c) 2012 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Efi Fogel <efifogel@gmail.com>
//             Utkarsh Khajuria <utkarshkhajuria55@gmail.com>

#ifndef CGAL_DRAW_ARRANGEMENT_2_H
#define CGAL_DRAW_ARRANGEMENT_2_H

#include <CGAL/Qt/Basic_viewer.h>

#include "CGAL/Bbox_2.h"

#ifdef DOXYGEN_RUNNING

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Draw
 *
 * The function opens a new window and draws `arr`, an object the type of which
 * is an instance of the `Arrangement_on_surface_2` class template. Parameters
 * of the drawing are taken from the graphics scene options `gso`.
 *
 * The arrangement is approximated within the visible part of the parameter
 * space, and it is approximated again whenever the camera changes, so that the
 * approximation error remains below the size of a pixel.
 *
 * A call to this function blocks the execution of the program until the drawing
 * window is closed. This function requires `CGAL_Qt6`, and is only available if
 * the macro `CGAL_USE_BASIC_VIEWER` is defined. Linking with the cmake target
 * `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition
 * `CGAL_USE_BASIC_VIEWER`.
 *
 * \tparam GeometryTraits a geometry traits type, a model of a 2D arrangement
 * geometry traits concept. It must support the approximation of points and
 * \f$x\f$-monotone curves; observe that not all geometry-traits models are
 * supported.
 * \tparam TopologyTraits a topology traits type, a model of the
 * `AosTopologyTraits` concept.
 * \tparam GSOptions a model of the `GraphicsSceneOptions` concept.
 *
 * \param arr the 2D arrangement to draw.
 * \param gso the graphics scene options.
 * \param title the title of the window.
 * \param bbox a bounding box in parameter space shown initially. If it is
 * empty (the default), a bounding box that contains the arrangement is used; if
 * the arrangement is induced by unbounded curves, this bounding box contains
 * the bounded parts of the curves. It is ignored if the arrangement is embedded
 * on the sphere, in which case the entire sphere is shown.
 *
 * \sa `AosTraits_2`
 * \sa `AosTopologyTraits`
 * \sa `GraphicsSceneOptions`
 */
template <typename GeometryTraits, typename TopologyTraits, typename GSOptions>
void draw(const Arrangement_on_surface_2<GeometryTraits, TopologyTraits>& arr,
          const GSOptions& gso,
          const char* title = "2D Arrangement on Surface Viewer",
          const Bbox_2& bbox = Bbox_2(0, 0, 0, 0));

/*! \ingroup PkgArrangementOnSurface2Draw
 *
 * A shortcut to `draw(arr, gso, title, bbox)`, where `gso` are graphics scene
 * options of type
 * `Graphics_scene_options<Aos, Aos::Vertex_const_handle, Aos::Halfedge_const_handle, Aos::Face_const_handle>`,
 * `Aos` is `Arrangement_on_surface_2<GeometryTraits, TopologyTraits>`,
 * vertices are red, edges are black, and faces are colored randomly.
 */
template <typename GeometryTraits, typename TopologyTraits>
void draw(const Arrangement_on_surface_2<GeometryTraits, TopologyTraits>& arr,
          const char* title = "2D Arrangement on Surface Viewer",
          const Bbox_2& bbox = Bbox_2(0, 0, 0, 0));

/*! \ingroup PkgArrangementOnSurface2Draw
 *
 * Adds the vertices, edges, and faces of `arr` to the graphics scene `gs`.
 * The scene is not cleared before it is filled; thus, several objects (e.g.,
 * several arrangements) can be added to the same scene and then drawn in a
 * single window by `draw_graphics_scene()`.
 *
 * The curves are approximated once, with an approximation error that is a
 * small fraction of the size of the bounding box, and the faces are
 * triangulated. (In contrast, `draw()` re-approximates the arrangement
 * whenever the camera changes.)
 *
 * \tparam GeometryTraits a geometry traits type, a model of a 2D arrangement
 * geometry traits concept. It must support the approximation of points and
 * \f$x\f$-monotone curves.
 * \tparam TopologyTraits a topology traits type, a model of the
 * `AosTopologyTraits` concept.
 * \tparam GSOptions a model of the `GraphicsSceneOptions` concept.
 *
 * \param arr the arrangement to add to the graphics scene.
 * \param gs the graphics scene to fill.
 * \param gso the graphics scene options.
 * \param bbox a bounding box in parameter space, to which the arrangement is
 * clipped. If it is empty (the default), a bounding box that contains the
 * arrangement is used; if the arrangement is induced by unbounded curves, this
 * bounding box contains the bounded parts of the curves. It is ignored if the
 * arrangement is embedded on the sphere, in which case the entire arrangement
 * is added.
 */
template <typename GeometryTraits, typename TopologyTraits, typename GSOptions>
void add_to_graphics_scene(const Arrangement_on_surface_2<GeometryTraits, TopologyTraits>& arr,
                           Graphics_scene& gs, const GSOptions& gso,
                           const Bbox_2& bbox = Bbox_2(0, 0, 0, 0));

/*! \ingroup PkgArrangementOnSurface2Draw
 *
 * A shortcut to `add_to_graphics_scene(arr, gs, gso, bbox)`, where `gso` are
 * graphics scene options of type
 * `Graphics_scene_options<Aos, Aos::Vertex_const_handle, Aos::Halfedge_const_handle, Aos::Face_const_handle>`,
 * `Aos` is `Arrangement_on_surface_2<GeometryTraits, TopologyTraits>`,
 * vertices are red, edges are black, and faces are colored randomly.
 */
template <typename GeometryTraits, typename TopologyTraits>
void add_to_graphics_scene(const Arrangement_on_surface_2<GeometryTraits, TopologyTraits>& arr, Graphics_scene& gs,
                           const Bbox_2& bbox = Bbox_2(0, 0, 0, 0));

} // namespace CGAL

#endif

#endif
