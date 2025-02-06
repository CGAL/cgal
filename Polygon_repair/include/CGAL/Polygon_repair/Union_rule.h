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

#ifndef CGAL_POLYGON_REPAIR_UNION_RULE_H
#define CGAL_POLYGON_REPAIR_UNION_RULE_H

#include <CGAL/license/Polygon_repair.h>

namespace CGAL {

namespace Polygon_repair {

/// \addtogroup PkgPolygonRepairRules
/// @{

/*!
  Tag class to select the %union rule when calling `CGAL::Polygon_repair::repair()`.
  The union rules are useful when given two or more similar valid polygons with holes.
The union rule results in areas that are contained in at least one of the input polygons with holes.
  */
  struct Union_rule {};

///@}

} // namespace Polygon_repair

} // namespace CGAL

#endif  // CGAL_POLYGON_REPAIR_UNION_RULE_H
