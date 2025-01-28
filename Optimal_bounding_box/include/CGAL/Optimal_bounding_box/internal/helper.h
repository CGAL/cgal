// Copyright (c) 2018-2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Konstantinos Katrioplas
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_OPTIMAL_BOUNDING_BOX_INTERNAL_HELPER_H
#define CGAL_OPTIMAL_BOUNDING_BOX_INTERNAL_HELPER_H

#include <CGAL/license/Optimal_bounding_box.h>

namespace CGAL {
namespace Optimal_bounding_box {
namespace internal {

template <typename Matrix>
Matrix transpose(const Matrix& m)
{
  Matrix tm;

  tm.set(0, 0, m(0, 0)); tm.set(0, 1, m(1, 0)); tm.set(0, 2, m(2, 0));
  tm.set(1, 0, m(0, 1)); tm.set(1, 1, m(1, 1)); tm.set(1, 2, m(2, 1));
  tm.set(2, 0, m(0, 2)); tm.set(2, 1, m(1, 2)); tm.set(2, 2, m(2, 2));

  return tm;
}

} // namespace internal
} // namespace Optimal_bounding_box
} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_BOX_INTERNAL_HELPER_H
