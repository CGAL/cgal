// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_INTERNAL_DOMAIN_CELL_TYPE_H
#define CGAL_ISOSURFACING_3_INTERNAL_DOMAIN_CELL_TYPE_H

#include <CGAL/license/Isosurfacing_3.h>

#include <limits>

namespace CGAL {
namespace Isosurfacing {

// Was supposed to check if an algorithm can handle a specific domain. Not used right now.
using Cell_type = std::size_t;

static constexpr Cell_type ANY_CELL = (std::numeric_limits<std::size_t>::max)();

static constexpr Cell_type POLYHEDRAL_CELL = (std::size_t(1) << 0);
static constexpr Cell_type TETRAHEDRAL_CELL = (std::size_t(1) << 1);
static constexpr Cell_type CUBICAL_CELL = (std::size_t(1) << 2);

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_DOMAIN_CELL_TYPE_H
