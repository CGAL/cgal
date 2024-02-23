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

#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef DOXYGEN_RUNNING

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2Draw
 *
 * opens a new window and draws `arr`, an instance of the `CGAL::Arrangement_2`
 * class template. A call to this function is blocking; that is, the program
 * continues only after the user closes the window. This function requires
 * `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is
 * defined.  Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link
 * with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.
 *
 * \tparam GeometryTraits_2 a geometry traits type, a model of a 2D arrangement
 * traits concept. At this point it must be an instance of either
 * `CGAL::Arr_segment_traits_2` or `CGAL::Arr_conic_traits_2`.
 *
 * \tparam Dcel the \dcel type, a model of the `ArrangementDcel` concept.
 *
 * \param arr the 2D arrangement to draw.
 * \param title the window title.
 *
 * \sa `ArrangementDcel`
 * \sa `ArrangementTraits_2`
 */
template <typename GeometryTraits_2, typename Dcel>
void draw(const Arrangement_2<GeometryTraits_2, Dcel>& arr,
          const char* title = "2D Arrangement Basic Viewer");

} /* namespace CGAL */

#endif

#endif
