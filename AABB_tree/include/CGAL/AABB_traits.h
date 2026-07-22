// Copyright (c) 2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Stéphane Tayeb, Pierre Alliez, Camille Wormser
//

#ifndef CGAL_AABB_TRAITS_H_
#define CGAL_AABB_TRAITS_H_

#define CGAL_DEPRECATED_HEADER "<CGAL/AABB_traits.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/AABB_traits_3.h>"
#include <CGAL/Installation/internal/deprecation_warning.h>

#ifndef CGAL_NO_DEPRECATED_CODE

#include <CGAL/license/AABB_tree.h>

#include <CGAL/AABB_traits_3.h>

/// \file AABB_traits.h


namespace CGAL
{



/// \addtogroup PkgAABBTreeRef
/// @{

/// template alias for backward compatibility
template<typename GeomTraits, typename AABBPrimitive, typename BboxMap = Default>
using AABB_traits = AABB_traits_3<GeomTraits, AABBPrimitive, BboxMap>;

///@}

} // namespace CGAL

#endif // CGAL_NO_DEPRECATED_CODE

#endif // CGAL_AABB_TRAITS_H_
