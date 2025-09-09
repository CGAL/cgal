// Copyright (c) 2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_WEAK_FIND_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_WEAK_FIND_H

#include <algorithm>
#include <iterator>
#include <memory>

namespace CGAL {
namespace STL_Extension {
namespace internal {

template<class Iter, class V>
auto weak_find(Iter begin, Iter end, const std::weak_ptr<V>& target)
{
  return std::find_if(begin, end, [&](const std::weak_ptr<V>& w) {
                                      auto sp = w.lock();
                                      auto tp = target.lock();
                                      return sp && tp && (sp == tp);
                                    });
}

} // namespace internal
} // namespace STL_Extension
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_WEAK_FIND_H */
