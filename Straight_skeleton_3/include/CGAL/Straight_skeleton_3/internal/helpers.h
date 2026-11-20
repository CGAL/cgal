// Copyright (c) 2026 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_HELPERS_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_HELPERS_H

#include <CGAL/license/Straight_skeleton_3.h>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {

// defining a sorter for maps and sets with T a pointer that sorts using ->id()
// Put higher Ids first
struct ID_sorter
{
  template <typename T>
  bool operator()(const T& a, const T& b) const
  {
    return a->id() > b->id();
  }
};

template <typename T>
using ID_set = std::set<T, ID_sorter>;

template <typename T>
using ID_map = std::map<T, T, ID_sorter>;

} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_HELPERS_H */
