// Copyright (c) 2024 GeometryFactory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_REPAIR_INTERSECTION_RULE_H
#define CGAL_POLYGON_REPAIR_INTERSECTION_RULE_H

#include <CGAL/license/Polygon_repair.h>

namespace CGAL {

namespace Polygon_repair {

/// \addtogroup PkgPolygonRepairRules
/// @{

/*!
  Tag class to select the %intersection rule when calling `CGAL::Polygon_repair::repair()`.
  The intersection rule are useful when given
two or more similar valid polygons with holes.
The intersection rule results in areas that are contained in all input polygons with holes.
  */
  struct Intersection_rule {};

///@}

} // namespace Polygon_repair

} // namespace CGAL

#endif  // CGAL_POLYGON_REPAIR_INTERSECTION_RULE_H
