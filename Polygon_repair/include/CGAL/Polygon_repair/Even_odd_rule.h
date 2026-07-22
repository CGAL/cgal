// Copyright (c) 2023 GeometryFactory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ken Arroyo Ohori

#ifndef CGAL_POLYGON_REPAIR_EVEN_ODD_RULE_H
#define CGAL_POLYGON_REPAIR_EVEN_ODD_RULE_H

#include <CGAL/license/Polygon_repair.h>

namespace CGAL {

namespace Polygon_repair {

/// \addtogroup PkgPolygonRepairRules
/// @{

/*!
  Tag class to select the even odd rule when calling `CGAL::Polygon_repair::repair()`.
  The even-odd rule results in areas that are alternately assigned as polygon
interiors and exterior/holes each time that an input edge is passed.
It does not distinguish between edges that are part of outer boundaries
from those of inner boundaries.
  */
  struct Even_odd_rule {};

///@}

} // namespace Polygon_repair

} // namespace CGAL

#endif  // CGAL_POLYGON_REPAIR_EVEN_ODD_RULE_H
