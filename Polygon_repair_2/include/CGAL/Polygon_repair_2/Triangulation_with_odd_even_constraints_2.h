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

#ifndef CGAL_TRIANGULATION_WITH_ODD_EVEN_CONSTRAINTS_2_H
#define CGAL_TRIANGULATION_WITH_ODD_EVEN_CONSTRAINTS_2_H

namespace CGAL {

/*! \ingroup PkgPolygonRepair2Ref
 *
 * The class `Triangulation_with_odd_even_constraints_2` builds on a constrained
 * triangulation to remove the parts of constraints that overlap an even number of times
 *
 * \tparam Triangulation_ must have support for constraints
 */
template <class Triangulation_>
class Triangulation_with_odd_even_constraints_2 {
public:
  /// \name Definition

  /// @{
  /// the triangulation class.
  typedef typename Triangulation_ Triangulation;

  /// handle to a vertex.
  typedef typename Triangulation_::Vertex_handle Vertex_handle;
  /// @}
  
  // iterator over interior faces.
  class Interior_faces_iterator : public Triangulation::All_faces_iterator {
    Interior_faces_iterator operator++();
    Interior_faces_iterator operator--();
  }

  // Add constraint from va to vb using the odd-even rule
  void odd_even_insert_constraint(Vertex_handle va, Vertex_handle vb);

  // Starts at an arbitrary interior face
  Interior_faces_iterator interior_faces_begin();

  // Past-the-end iterator
  Interior_faces_iterator interior_faces_end();
}

} // namespace CGAL

#endif // CGAL_TRIANGULATION_WITH_ODD_EVEN_CONSTRAINTS_2_H