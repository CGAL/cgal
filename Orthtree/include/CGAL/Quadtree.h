// Copyright (c) 2020  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_QUADTREE_H
#define CGAL_QUADTREE_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Orthtree.h>
#include <CGAL/Orthtree_traits_2.h>

namespace CGAL {

template <typename GeomTraits, typename PointRange,
          typename PointMap = Identity_property_map
         <typename std::iterator_traits<typename PointRange::iterator>::value_type> >
using Quadtree = Orthtree<Orthtree_traits_2<GeomTraits>, PointRange, PointMap>;

} // namespace CGAL


#endif // CGAL_OCTREE_H
