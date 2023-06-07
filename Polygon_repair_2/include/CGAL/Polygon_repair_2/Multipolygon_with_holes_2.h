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

#ifndef CGAL_MULTIPOLYGON_WITH_HOLES_2_H
#define CGAL_MULTIPOLYGON_WITH_HOLES_2_H

#include <CGAL/Polygon_with_holes_2.h>

namespace CGAL {

/*! \ingroup PkgPolygonRepair2Ref
 *
 * The class `Multipolygon_with_holes_2` models the concept
 * `MultiPolygonWithHoles_2`.
 * It is parameterized with a type `Polygon` used to define the exposed
 * type `Polygon_with_holes_2`. This type represents each polygon.
 *
 * \tparam Polygon_ must have input and output operators.
 *
 * \cgalModels `MultipolygonWithHoles_2`
 */
template <class PolygonTraits_,
          class PolygonContainer_ = std::vector<typename Kernel::Point_2>>
class Multipolygon_with_holes_2 {
public:
  /// \name Definition

  /// polygon with holes type
  typedef CGAL::Polygon_with_holes_2<PolygonTraits_, PolygonContainer_> Polygon_with_holes_2;
  
  typedef std::deque<Polygon_with_holes_2> Polygons_container;

  typedef typename Polygons_container::iterator Polygon_iterator;
  typedef typename Polygons_container::const_iterator Polygon_const_iterator;

  typedef unsigned int Size;
  

  Multipolygon_with_holes_2() {}

  template <typename PolygonsInputIterator>
  Multipolygon_with_holes_2(PolygonsInputIterator p_begin,
                            PolygonsInputIterator p_end) :
    m_polygons(p_begin, p_end)
  {}

  Polygons_container& polygons() { return m_polygons; }

  const Polygons_container& polygons() const { return m_polygons; }

  Polygon_iterator polygons_begin() { return m_polygons.begin(); }

  Polygon_iterator polygons_end() { return m_polygons.end(); }

  Polygon_const_iterator polygons_begin() const { return m_polygons.begin(); }

  Polygon_const_iterator polygons_end() const { return m_polygons.end(); }

  void add_polygon(const Polygon_with_holes_2& pgn) { m_polygons.push_back(pgn); }

  void add_polygon(Polygon_with_holes_2&& pgn) { m_polygons.emplace_back(std::move(pgn)); }

  void erase_polygon(Polygon_iterator pit) { m_polygons.erase(pit); }

  void clear() { m_polygons.clear(); }

  Size number_of_polygons() const { return static_cast<Size>(m_polygons.size()); }

protected:
  Polygons_container m_polygons;
};


} //namespace CGAL

#endif // CGAL_MULTIPOLYGON_WITH_HOLES_2_H
