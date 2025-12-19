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
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_AABB_TRIANGLE_PRIMITIVE_H_
#define CGAL_AABB_TRIANGLE_PRIMITIVE_H_

#include <CGAL/license/AABB_tree.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/AABB_triangle_primitive.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/AABB_triangle_primitive_3.h>"
#include <CGAL/Installation/internal/deprecation_warning.h>

#ifndef CGAL_NO_DEPRECATED_CODE

#include <CGAL/AABB_triangle_primitive_3.h>

/// \file AABB_triangle_primitive.h

namespace CGAL {


/// \addtogroup PkgAABBTreeRef
/// @{

/// template alias for backward compatibility

template < class GeomTraits,
           class Iterator,
           class CacheDatum=Tag_false>
using  AABB_triangle_primitive = AABB_triangle_primitive_3<GeomTraits, Iterator, CacheDatum>;

///@}

} //CGAL namespace

#endif  // CGAL_NO_DEPRECATED_CODE

#endif //CGAL_AABB_TRIANGLE_PRIMITIVE_H_
