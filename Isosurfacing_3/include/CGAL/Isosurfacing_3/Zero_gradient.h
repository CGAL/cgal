// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_ZERO_GRADIENT_H
#define CGAL_ISOSURFACING_3_ZERO_GRADIENT_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Origin.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domain_helpers_grp
 *
 * \brief Class template for a gradient that equals zero everywhere.
 *
 * \tparam P the point type
 *
 * \details This gradient function can be useful when using Marching Cubes, which does not require the domain to have a gradient.
 */
struct Zero_gradient
{
  /**
   * \return the null vector
   */
  template <typename P>
  CGAL::Null_vector operator()(const P&) const
  {
    return CGAL::NULL_VECTOR;
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_ZERO_GRADIENT_H
