// Copyright (c) 2012
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
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

 * The function opens a new window and draws `arr`, an object the type of which
 * is an instance of the `Arrangement_on_surface_2` class template. Parameters
 * of the drawing are taken from the optional graphics scene options parameter.
 *
 * A call to this function blocks the execution of the program until the drawing
 * window is closed. This function requires `CGAL_Qt6`, and is only available if
 * the macro `CGAL_USE_BASIC_VIEWER` is defined. Linking with the cmake target
 * `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition
 * `CGAL_USE_BASIC_VIEWER`.
 *
 * \tparam GeometryTraits a geometry traits type, a model of a 2D arrangement
 * geometry traits concept. Observe that not all geometery-traits models are
 * supported.
 *
 * \tparam TopologyTraits a topology traits type, a model of the
 * `AosTopologyTraits` concept.
 *
 * \tparam GSOptions a model of `GraphicsSceneOptions` concept.
 *
 * \param arr the 2D arrangement to draw.
 * \param bbox a bounding box in parameter space.
 * \param gso the graphics scene options.
 * \param title the optional title of the window.
 *
 * \sa `AosTraits_2`
 * \sa `AosTopologyTraits`
 * \sa `GraphicsSceneOptions`
 */

template <typename GeometryTraits, typename TopologyTraits>
void draw(const Arrangement_on_surface_2<GeometryTraits, TopologyTraits>& arr,
          const Bbox_2& bbox, const GSOptions& gso,
          const char* title = "2D Arrangement on Surface");

/*! \ingroup PkgArrangementOnSurface2Draw
 *
 * A shortcut to
 * `draw(arr, bbox, Graphics_scene_options<Aos, Aos::Vertex_const_handle, Aos::Halfedge_const_handle, Aos::Face_const_handle>{})`,
 * where `Aos` is `Arrangement_on_surface_2<GeometryTraits, TopologyTraits>`.
 */

template <typename GeometryTraits, typename TopologyTraits>
void draw(const Arrangement_on_surface_2<GeometryTraits, TopologyTraits>& arr,
          const Bbox_2& bbox, const char* title = "2D Arrangement on Surface");

/*! \ingroup PkgArrangementOnSurface2Draw
 *
 * Similar to `draw(arr, bbox, gso)`, where the bounding box `bbox` is
 * computed to bound all points and curves of the arrangement in parameter
 * space.
 */

template <typename GeometryTraits, typename TopologyTraits>
void draw(const Arrangement_on_surface_2<GeometryTraits, TopologyTraits>& arr,
          const GSOptions& gso, const char* title = "2D Arrangement on Surface");

/*! \ingroup PkgArrangementOnSurface2Draw
 *
 * A shortcut to
 * `draw(arr, Graphics_scene_options<Aos, Aos::Vertex_const_handle, Aos::Halfedge_const_handle, Aos::Face_const_handle>{})`,
 * where `Aos` is `Arrangement_on_surface_2<GeometryTraits, TopologyTraits>`.
 */

template <typename GeometryTraits, typename TopologyTraits>
void draw(const Arrangement_on_surface_2<GeometryTraits, TopologyTraits>& arr,
          const char* title = "2D Arrangement on Surface");

/*! \ingroup PkgArrangementOnSurface2Draw
 *
 * Adds the vertices, edges, and faces of `arr` to the graphics scene `gs`.
 * The scene is not cleared before it is filled.
 *
 * \tparam GeometryTraits a geometry traits type, a model of a 2D arrangement
 * geometry traits concept.
 * \tparam TopologyTraits a topology traits type, a model of the
 * `AosTopologyTraits` concept.
 * \tparam GSOptions a model of the `GraphicsSceneOptions` concept.
 *
 * \param arr the arrangement to add to the graphics scene.
 * \param gs the graphics scene to fill.
 * \param gso the graphics scene options.
 */
template <typename GeometryTraits, typename TopologyTraits, typename GSOptions>
void add_to_graphics_scene(const Arrangement_on_surface_2<GeometryTraits, TopologyTraits>& arr,
                           Graphics_scene& gs, const GSOptions& gso);

/*! \ingroup PkgArrangementOnSurface2Draw
 *
 * A shortcut to
 * `add_to_graphics_scene(arr, gs, Graphics_scene_options<Aos, Aos::Vertex_const_handle, Aos::Halfedge_const_handle, Aos::Face_const_handle>{})`,
 * where `Aos` is `Arrangement_on_surface_2<GeometryTraits, TopologyTraits>`.
 *
 * \param arr the arrangement to add to the graphics scene.
 * \param gs the graphics scene to fill.
 */
template <typename GeometryTraits, typename TopologyTraits>
void add_to_graphics_scene(const Arrangement_on_surface_2<GeometryTraits, TopologyTraits>& arr, Graphics_scene& gs);

} // namespace CGAL

#endif

#endif
