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
template <class PolygonTraits_, class PolygonContainer_>
Multipolygon_with_holes_2<class PolygonTraits_, class PolygonContainer_> repair(const Polygon_2<PolygonTraits_, PolygonContainer_>& p);

/// \ingroup PkgPolygonRepair2Functions
/// Repair a polygon with holes
template <class PolygonTraits_, class PolygonContainer_>
Multipolygon_with_holes_2<class PolygonTraits_, class PolygonContainer_> repair(const Polygon_with_holes_2<PolygonTraits_, PolygonContainer_>& p);

/// \ingroup PkgPolygonRepair2Functions
/// Repair a multipolygon with holes
Multipolygon_with_holes_2<class PolygonTraits_, class PolygonContainer_> repair(const Multipolygon_with_holes_2& p);

/*! \ingroup PkgPolygonRepair2Ref
 *
 * The class `Polygon_repair_2` builds on a constrained
 * triangulation to remove the parts of constraints that overlap an even number of times
 */
template <class PolygonTraits_, class PolygonContainer_>
class Polygon_repair_2 {
public:

  /// \name Creation
  Polygon_repair_2();

	/// \name Modifiers
  /// @{

	// Add edges of the polygon to the triangulation
	void add_to_triangulation(const Polygon_2<PolygonTraits_, PolygonContainer_>& p);

	// Add edges of the polygon to the triangulation
	void add_to_triangulation(const Polygon_with_holes_2<PolygonTraits_, PolygonContainer_>& p);

	// Add edges of the polygon to the triangulation
	void add_to_triangulation(const Multipolygon_with_holes_2<PolygonTraits_, PolygonContainer_>& p);

	// Label triangles in triangulation as inside or outside the polygon
	void label_triangulation();

  // Reconstruct multipolygon based on the triangles labelled as inside the polygon
  void reconstruct_polygon();

  // Erases the triangulation.
  void clear();

	/// @}

	/// \name Access Functions
  /// @{

  Triangulation_with_odd_even_constraints_2 triangulation();

  Multipolygon_with_holes_2<PolygonTraits_, PolygonContainer_> multipolygon();

  /// @}
	

protected:
	Triangulation_with_odd_even_constraints_2 t;
  Multipolygon_with_holes_2<PolygonTraits_, PolygonContainer_> mp;
}

} // namespace Polygon_repair_2
} // namespace CGAL

#endif // CGAL_POLYGON_REPAIR_2_H
