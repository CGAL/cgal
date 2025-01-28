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

/// \addtogroup PkgPolygonRepairRef
/// @{

/*!
  Tag class to select the even odd rule when calling `CGAL::Polygon_repair::repair()`.
  */
  struct Even_odd_rule {};

///@}

} // namespace Polygon_repair

} // namespace CGAL

#endif  // CGAL_POLYGON_REPAIR_EVEN_ODD_RULE_H
