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
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Polygon_repair_2/Multipolygon_with_holes_2.h>
#include <CGAL/Polygon_repair_2/Triangulation_with_odd_even_constraints_2.h>

namespace CGAL {

namespace Polygon_repair_2 {

template <class Kernel, class PolygonContainer>
class Polygon_repair_2;

/// \ingroup PkgPolygonRepair2Functions
/// Repair a polygon without holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair(const CGAL::Polygon_2<Kernel, PolygonContainer>& p) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation(p);
  pr.label_triangulation();
  pr.reconstruct_polygon();
  return pr.multipolygon();
}

/// \ingroup PkgPolygonRepair2Functions
/// Repair a polygon with holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair(const CGAL::Polygon_with_holes_2<Kernel, PolygonContainer>& p) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation(p);
  pr.label_triangulation();
  pr.reconstruct_polygon();
  return pr.multipolygon();
}

/// \ingroup PkgPolygonRepair2Functions
/// Repair a multipolygon with holes
template <class Kernel, class PolygonContainer>
Multipolygon_with_holes_2<Kernel, PolygonContainer> repair(const Multipolygon_with_holes_2<Kernel, PolygonContainer>& mp) {
  CGAL::Polygon_repair_2::Polygon_repair_2<Kernel, PolygonContainer> pr;
  pr.add_to_triangulation(mp);
  pr.label_triangulation();
  pr.reconstruct_polygon();
  return pr.multipolygon();
}

/*! \ingroup PkgPolygonRepair2Ref
 *
 * The class `Polygon_repair_2` builds on a constrained
 * triangulation to remove the parts of constraints that overlap an even number of times
 */
template <class Kernel = CGAL::Exact_predicates_inexact_constructions_kernel,
          class PolygonContainer = std::vector<typename Kernel::Point_2>>
class Polygon_repair_2 {
public:

  struct Repair_face_info {
    bool processed;
    bool interior;
    Repair_face_info() {
      processed = false;
      interior = false;
    }
  };
  typedef CGAL::Triangulation_vertex_base_2<Kernel> Vertex_base;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Face_base;
  typedef CGAL::Triangulation_face_base_with_info_2<Repair_face_info, Kernel, Face_base> Face_base_with_repair_info;
  typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base_with_repair_info> Triangulation_data_structure;
  typedef CGAL::Exact_predicates_tag Tag; // assumed for now
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Triangulation_data_structure, Tag> Constrained_Delaunay_triangulation;
  typedef Triangulation_with_odd_even_constraints_2<Constrained_Delaunay_triangulation> Triangulation;

  /// \name Creation
  Polygon_repair_2() {

  }

	/// \name Modifiers
  /// @{

	// Add edges of the polygon to the triangulation
	void add_to_triangulation(const Polygon_2<Kernel, PolygonContainer>& p) {
    for (auto const& e: p.edges()) {
      t.odd_even_insert_constraint(e.source(), e.target());
    }
  }

	// Add edges of the polygon to the triangulation
	void add_to_triangulation(const Polygon_with_holes_2<Kernel, PolygonContainer>& p) {
    add_to_triangulation(p.outer_boundary());
    for (auto const& h: p.holes()) {
       add_to_triangulation(h);
    }
  }

	// Add edges of the polygon to the triangulation
	void add_to_triangulation(const Multipolygon_with_holes_2<Kernel, PolygonContainer>& mp) {
    for (auto const &p: mp.polygons()) {
      add_to_triangulation(p);
    }
  }

	// Label triangles in triangulation as inside or outside the polygon
	void label_triangulation() {

  }

  // Reconstruct multipolygon based on the triangles labelled as inside the polygon
  void reconstruct_polygon() {
    
  }

  // Erases the triangulation.
  void clear() {
    t.clear();
  }

	/// @}

	/// \name Access Functions
  /// @{

  Triangulation_with_odd_even_constraints_2<Triangulation>& triangulation() {
    return t;
  }

  Multipolygon_with_holes_2<Kernel, PolygonContainer> multipolygon() {

  }

  /// @}
	

protected:
	Triangulation_with_odd_even_constraints_2<Triangulation> t;
  Multipolygon_with_holes_2<Kernel, PolygonContainer> mp;
};

} // namespace Polygon_repair_2
} // namespace CGAL

#endif // CGAL_POLYGON_REPAIR_2_H
