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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_DOMAIN_CELL_TYPE
#define CGAL_ISOSURFACING_3_INTERNAL_DOMAIN_CELL_TYPE

#include <CGAL/license/Isosurfacing_3.h>

#include <limits>

namespace CGAL {
namespace Isosurfacing {

// Was supposed to check if an algorithm can handle a specific domain. Not used right now. @fixme remove
using Cell_type = std::size_t;

static constexpr Cell_type ANY_CELL = (std::numeric_limits<Cell_type>::max)();

static constexpr Cell_type POLYHERDAL_CELL = (Cell_type(1) << 0);
static constexpr Cell_type TETRAHEDRAL_CELL = (Cell_type(1) << 1);
static constexpr Cell_type CUBICAL_CELL = (Cell_type(1) << 2);

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_DOMAIN_CELL_TYPE
