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

#ifndef CGAL_DRAW_ARRANGEMENT_2_H
#define CGAL_DRAW_ARRANGEMENT_2_H

#include <CGAL/Qt/Basic_viewer.h>

#ifdef DOXYGEN_RUNNING

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Draw

  opens a new window and draws `arr`, an instance of the `CGAL::Arrangement_2` class template. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.


\tparam GeometryTraits_2 a geometry traits type, a model of a 2D arrangement traits concept. At this point it must be an instance of either `CGAL::Arr_segment_traits_2` or `CGAL::Arr_conic_traits_2`.
\tparam Dcel the \dcel type, a model of the `ArrangementDcel` concept.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param arr the 2D arrangement to draw.
\param gso the graphics scene options parameter.

\sa `ArrangementDcel`
\sa `ArrangementTraits_2`
*/
template <typename GeometryTraits_2, typename Dcel, typename GSOptions>
void draw(const Arrangement_2<GeometryTraits_2, Dcel>& arr, const GSOptions& gso);

/*! \ingroup PkgArrangementOnSurface2Draw

  A shortcut to `CGAL::draw(arr, Graphics_scene_options{})`.
*/
template <typename GeometryTraits_2, typename Dcel>
void draw(const Arrangement_2<GeometryTraits_2, Dcel>& arr);

/*! \ingroup PkgArrangementOnSurface2Draw

adds the vertices, edges and faces of `arr` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam GeometryTraits_2 a geometry traits type, a model of a 2D arrangement traits concept. At this point it must be an instance of either `CGAL::Arr_segment_traits_2` or `CGAL::Arr_conic_traits_2`.
\tparam Dcel the \dcel type, a model of the `ArrangementDcel` concept.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param arr the 2D arrangement to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.
 */
template <typename GeometryTraits_2, typename Dcel, typename GSOptions>
void add_to_graphics_scene(const Arrangement_2<GeometryTraits_2, Dcel>& arr,
                           CGAL::Graphics_scene& gs, const GSOptions& gso);

/*! \ingroup PkgArrangementOnSurface2Draw

  A shortcut to `CGAL::add_to_graphics_scene(arr, gs, Graphics_scene_options{})`.
 */
template <typename GeometryTraits_2, typename Dcel>
void add_to_graphics_scene(const Arrangement_2<GeometryTraits_2, Dcel>& arr,
                           CGAL::Graphics_scene& gs);

} /* namespace CGAL */

#endif

#endif
