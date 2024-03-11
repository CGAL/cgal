// Copyright (c) 2005
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ken Arroyo Ohori <k.ohori@tudelft.nl>

#ifndef CGAL_MULTIPOLYGON_WITH_HOLES_2_H
#define CGAL_MULTIPOLYGON_WITH_HOLES_2_H

#include <CGAL/Polygon_with_holes_2.h>

namespace CGAL {

/*! \ingroup PkgPolygon2Ref
 *
 * The class `Multipolygon_with_holes_2` models the concept `MultipolygonWithHoles_2`.
 * It is parameterized with two types (`Kernel` and `Container`) that are used to instantiate
 * the types `Polygon_2<Kernel,Container>` and `Polygon_with_holes_2<Kernel,Container>`.
 * The latter is used to represent each polygon with holes. The former is converted to the latter.
 *
 * \cgalModels{MultipolygonWithHoles_2}
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

  using Polygon_with_holes_container = std::deque<Polygon_with_holes_2>;

  using Polygon_with_holes_iterator = typename Polygon_with_holes_container::iterator;
  using Polygon_with_holes_const_iterator = typename Polygon_with_holes_container::const_iterator;

  /// the size type
  using Size = unsigned int;

  /*! %Default constructor. */
  Multipolygon_with_holes_2() {}

  /*! Constructor from polygons. */
  template <typename PolygonsInputIterator>
  Multipolygon_with_holes_2(PolygonsInputIterator p_begin,
                            PolygonsInputIterator p_end) :
    m_polygons(p_begin, p_end)
  {}

  Polygon_with_holes_container& polygons_with_holes() { return m_polygons; }

  const Polygon_with_holes_container& polygons_with_holes() const { return m_polygons; }

  Polygon_with_holes_iterator polygons_with_holes_begin() { return m_polygons.begin(); }

  Polygon_with_holes_iterator polygons_with_holes_end() { return m_polygons.end(); }

  Polygon_with_holes_const_iterator polygons_with_holes_begin() const { return m_polygons.begin(); }

  Polygon_with_holes_const_iterator polygons_with_holes_end() const { return m_polygons.end(); }

  void add_polygon(const Polygon_2& pgn) { m_polygons.push_back(Polygon_with_holes_2(pgn)); }

  void add_polygon_with_holes(const Polygon_with_holes_2& pgn) { m_polygons.push_back(pgn); }

  void add_polygon_with_holes(Polygon_with_holes_2&& pgn) { m_polygons.emplace_back(std::move(pgn)); }

  void erase_polygon_with_holes(Polygon_with_holes_iterator pit) { m_polygons.erase(pit); }

  void clear() { m_polygons.clear(); }

  Size number_of_polygons_with_holes() const { return static_cast<Size>(m_polygons.size()); }

protected:
  Polygon_with_holes_container m_polygons;
};


template <class Kernel_, class Container_>
bool operator==(const Multipolygon_with_holes_2<Kernel_, Container_>& p1,
                const Multipolygon_with_holes_2<Kernel_, Container_>& p2)
{
  typedef typename
    Multipolygon_with_holes_2<Kernel_, Container_>::Polygon_with_holes_const_iterator HCI;
  typedef CGAL::Polygon_with_holes_2<Kernel_, Container_> Polygon_2;
  if(&p1 == &p2)
    return (true);

  if(p1.number_of_polygons_with_holes() != p2.number_of_polygons_with_holes())
    return (false);

  std::list<Polygon_2> tmp_list(p2.polygons_with_holes_begin(), p2.polygons_with_holes_end());

  HCI i = p1.polygons_with_holes_begin();
  for(; i!= p1.polygons_with_holes_end(); ++i)
  {
    typename std::list<Polygon_2>::iterator j =
      (std::find(tmp_list.begin(), tmp_list.end(), *i));

    if(j == tmp_list.end())
      return (false);

    tmp_list.erase(j);
  }


  CGAL_assertion(tmp_list.empty());
  return (true);
}

template <class Kernel_, class Container_>
inline bool operator!=(const Multipolygon_with_holes_2<Kernel_, Container_>& p1,
                       const Multipolygon_with_holes_2<Kernel_, Container_>& p2)
{
  return (!(p1==p2));
}
/*!
inserts a multipolygon with holes to the output stream `os`.

An \ascii and a binary format exist. The format can be selected with
the \cgal modifiers for streams, `set_ascii_mode()` and `set_binary_mode()`,
respectively. The modifier `set_pretty_mode()` can be used to allow for (a
few) structuring comments in the output. Otherwise, the output would
be free of comments. The default for writing is \ascii without comments.

The number of polygons is written followed by the polygons. For each polygon,
the number of points of the outer boundary is written followed by the
points themselves in counterclockwise order. Then, the number of holes
is written, and for each hole, the number of points on its outer
boundary is written followed by the points themselves in clockwise
order.

\relates Multipolygon_with_holes_2
*/
template <class Kernel, class Container>
std::ostream& operator<<(std::ostream& os,
                         const Multipolygon_with_holes_2<Kernel, Container>& mp) {
  typename Multipolygon_with_holes_2<Kernel, Container>::Polygon_with_holes_const_iterator i;

  switch(IO::get_mode(os)) {
    case IO::ASCII :
      os << mp.number_of_polygons_with_holes() << ' ';
      for (i = mp.polygons_with_holes_begin(); i != mp.polygons_with_holes_end(); ++i) {
        os << *i << ' ';
      }
      return os;

    case IO::BINARY :
      os << mp.number_of_polygons_with_holes();
      for (i = mp.polygons_with_holes_begin(); i != mp.polygons_with_holes_end(); ++i) {
        os << *i ;
      }
      return os;

    default:
      os << "Multipolygon_with_holes_2(" << std::endl;
      for (i = mp.polygons_with_holes_begin(); i != mp.polygons_with_holes_end(); ++i) {
        os << " " << *i << std::endl;
      }

      os << ")" << std::endl;
      return os;
  }
}


} //namespace CGAL

#endif // CGAL_MULTIPOLYGON_WITH_HOLES_2_H
