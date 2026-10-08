// Copyright (c) 2025 Geometry Factory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Léo Valque

#ifndef CGAL_SNAP_ROUNDING_2_TRAITS_H
#define CGAL_SNAP_ROUNDING_2_TRAITS_H

#include <CGAL/license/Snap_rounding_2.h>

#include <CGAL/Hot_pixel_snap_rounding_traits_2.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Snap_rounding_traits_2.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Hot_pixel_snap_rounding_traits_2.h>"
#include <CGAL/Installation/internal/deprecation_warning.h>

namespace CGAL {

#ifndef CGAL_NO_DEPRECATED_CODE

/// \deprecated Use Hot_pixel_snap_rounding_traits_2 instead
template<class K>
using Snap_rounding_traits_2 = Hot_pixel_snap_rounding_traits_2<K>;
#endif // CGAL_NO_DEPRECATED_CODE

} //namespace CGAL

#endif // CGAL_ISR_2_TRAITS_H
