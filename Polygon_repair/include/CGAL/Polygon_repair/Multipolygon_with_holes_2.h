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

#ifndef CGAL_MULTIPOLYGON_WITH_HOLES_2_H
#define CGAL_MULTIPOLYGON_WITH_HOLES_2_H

#include <CGAL/license/Polygon_repair.h>

#include <CGAL/Polygon_with_holes_2.h>

namespace CGAL {

/*! \ingroup PkgPolygonRepairRef
 *
 * The class `Multipolygon_with_holes_2` models the concept `MultipolygonWithHoles_2`.
 * It is parameterized with two types (`Kernel` and `Container`) that are used to instantiate
 * the types `Polygon_2<Kernel,Container>` and `Polygon_with_holes_2<Kernel,Container>`.
 * The latter is used to represent each polygon with holes. The former is converted to the latter.
 *
 * \cgalModels `MultipolygonWithHoles_2`
 */
template <class Kernel,
          class Container = std::vector<typename Kernel::Point_2>>
class Multipolygon_with_holes_2 {
public:
  /// \name Definition

  /// @{

  /// polygon type
  using Polygon_2 = CGAL::Polygon_2<Kernel, Container>;

  /// polygon with holes type
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel, Container>;

  /// @}

  using Polygons_container = std::deque<Polygon_with_holes_2>;

  using Polygon_iterator = typename Polygons_container::iterator;
  using Polygon_const_iterator = typename Polygons_container::const_iterator;

  using Size = unsigned int;

  /*! %Default constructor. */
  Multipolygon_with_holes_2() {}

  /*! Constructor from polygons. */
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

  void add_polygon(const Polygon_2& pgn) { m_polygons.push_back(Polygon_with_holes_2(pgn)); }

  void add_polygon(const Polygon_with_holes_2& pgn) { m_polygons.push_back(pgn); }

  void add_polygon(Polygon_with_holes_2&& pgn) { m_polygons.emplace_back(std::move(pgn)); }

  void erase_polygon(Polygon_iterator pit) { m_polygons.erase(pit); }

  void clear() { m_polygons.clear(); }

  Size number_of_polygons() const { return static_cast<Size>(m_polygons.size()); }

protected:
  Polygons_container m_polygons;
};

/*!
exports a multipolygon with holes to the output stream `os`.

An \ascii and a binary format exist. The format can be selected with
the \cgal modifiers for streams, `set_ascii_mode()` and `set_binary_mode()`,
respectively. The modifier `set_pretty_mode()` can be used to allow for (a
few) structuring comments in the output. Otherwise, the output would
be free of comments. The default for writing is \ascii without comments.

The number of polygons is exported followed by the polygons. For each polygon,
the number of points of the outer boundary is exported followed by the
points themselves in counterclockwise order. Then, the number of holes
is exported, and for each hole, the number of points on its outer
boundary is exported followed by the points themselves in clockwise
order.

\relates Multipolygon_with_holes_2
*/
template <class Kernel, class Container>
std::ostream& operator<<(std::ostream& os,
                         const Multipolygon_with_holes_2<Kernel, Container>& mp) {
  typename Multipolygon_with_holes_2<Kernel, Container>::Polygon_const_iterator i;

  switch(IO::get_mode(os)) {
    case IO::ASCII :
      os << mp.number_of_polygons() << ' ';
      for (i = mp.polygons_begin(); i != mp.polygons_end(); ++i) {
        os << *i << ' ';
      }
      return os;

    case IO::BINARY :
      os << mp.number_of_polygons();
      for (i = mp.polygons_begin(); i != mp.polygons_end(); ++i) {
        os << *i ;
      }
      return os;

    default:
      os << "Multipolygon_with_holes_2(" << std::endl;
      for (i = mp.polygons_begin(); i != mp.polygons_end(); ++i) {
        os << " " << *i << std::endl;
      }

      os << ")" << std::endl;
      return os;
  }
}


} //namespace CGAL

#endif // CGAL_MULTIPOLYGON_WITH_HOLES_2_H
