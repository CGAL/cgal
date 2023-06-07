// Copyright (c) 2023 GeometryFactory. All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ken Arroyo Ohori

#ifndef CGAL_POLYGON_REPAIR_2_H
#define CGAL_POLYGON_REPAIR_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_repair_2/Multipolygon_with_holes_2.h>
#include <CGAL/Polygon_repair_2/Triangulation_with_odd_even_constraints_2.h>

namespace CGAL {

namespace Polygon_repair_2 {

/// \ingroup PkgPolygonRepair2Functions
/// Repair a polygon without holes
Multipolygon_with_holes_2 repair(const Polygon_2& p);

/// Repair a polygon with holes
Multipolygon_with_holes_2 repair(const Polygon_with_holes_2& p);

/// Repair a polygon with holes
Multipolygon_with_holes_2 repair(const Multipolygon_with_holes_2& p);

class Polygon_repair_2 {
  public:

  	/// \name Modifiers
	/// @{

  	// Add edges of the polygon to the triangulation
  	void add_to_triangulation(const Polygon_2& p);

  	// Add edges of the polygon to the triangulation
  	void add_to_triangulation(const Polygon_with_holes_2& p);

  	// Add edges of the polygon to the triangulation
  	void add_to_triangulation(const Multipolygon_with_holes_2& p);

  	// Label triangles in triangulation as inside or outside the polygon
  	void label_triangulation();

  	/// @}

  	// Reconstruct multipolygon based on the triangles labelled as inside the polygon
  	Multipolygon_with_holes_2 reconstruct_polygon();

  	// Erases the triangulation.
  	void clear();

  protected:
  	Triangulation_with_odd_even_constraints_2 triangulation;
}

} // namespace Polygon_repair_2
} // namespace CGAL

#endif // CGAL_POLYGON_REPAIR_2_H
