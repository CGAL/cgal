// Copyright (c) 2005
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GENERAL_POLYGON_WITH_HOLES_2_H
#define CGAL_GENERAL_POLYGON_WITH_HOLES_2_H

#include <list>
#include <iostream>
#include <CGAL/IO/io.h>

namespace CGAL {

/*!
\ingroup PkgPolygon2

The class `General_polygon_with_holes_2` models the concept
`GeneralPolygonWithHoles_2`. It represents a general polygon with
holes. It is parameterized with a type `Polygon` used to define
the exposed type `General_polygon_2`. This type represents the
outer boundary of the general polygon and the outer boundaries of
each hole.

\cgalModels `GeneralPolygonWithHoles_2`

*/
template <class Polygon_>
class General_polygon_with_holes_2
{
public:

/// \name Definition

/// @{
  typedef Polygon_							General_polygon_2;
/// @}

  typedef std::list<Polygon_>                         Holes_container;

  typedef typename Holes_container::iterator          Hole_iterator;
  typedef typename Holes_container::const_iterator    Hole_const_iterator;

  typedef unsigned int                                 Size;

  General_polygon_with_holes_2() : m_pgn()
  {}


  explicit General_polygon_with_holes_2(const General_polygon_2& pgn_boundary)
  : m_pgn(pgn_boundary)
  {}


  template <class HolesInputIterator>
  General_polygon_with_holes_2(const General_polygon_2& pgn_boundary,
                       HolesInputIterator h_begin,
                       HolesInputIterator h_end) : m_pgn(pgn_boundary),
                                                   m_holes(h_begin, h_end)
  {}

  Hole_iterator holes_begin()
  {
    return m_holes.begin();
  }

  Hole_iterator holes_end()
  {
    return m_holes.end();
  }

  Hole_const_iterator holes_begin() const
  {
    return m_holes.begin();
  }

  Hole_const_iterator holes_end() const
  {
    return m_holes.end();
  }

  bool is_unbounded() const
  {
    return m_pgn.is_empty();
  }

  General_polygon_2& outer_boundary()
  {
    return m_pgn;
  }

  const General_polygon_2& outer_boundary() const
  {
    return m_pgn;
  }

  void add_hole(const General_polygon_2& pgn_hole)
  {
    m_holes.push_back(pgn_hole);
  }

  void erase_hole(Hole_iterator hit)
  {
    m_holes.erase(hit);
  }

  bool has_holes() const
  {
    return (!m_holes.empty());
  }

  Size number_of_holes() const
  {
    return static_cast<Size>(m_holes.size());
  }

  void clear()
  {
    m_pgn.clear();
    m_holes.clear();
  }

  bool is_plane() const
  {
    return (m_pgn.is_empty() && m_holes.empty());
  }



protected:

  General_polygon_2           m_pgn;
  Holes_container            m_holes;
};


//-----------------------------------------------------------------------//
//                          operator<<
//-----------------------------------------------------------------------//
/*!
This operator exports a General_polygon_with_holes_2 to the output stream `out`.

An ASCII and a binary format exist. The format can be selected with
the \cgal modifiers for streams, `set_ascii_mode(0` and `set_binary_mode()`
respectively. The modifier `set_pretty_mode()` can be used to allow for (a
few) structuring comments in the output. Otherwise, the output would
be free of comments. The default for writing is ASCII without
comments.

The number of curves of the outer boundary is exported followed by the
curves themselves in counterclockwise order. Then, the number of holes
is exported, and for each hole, the number of curves on its outer
boundary is exported followed by the curves themselves in clockwise
order.

\relates General_polygon_with_holes_2
*/
template <class Polygon_>
std::ostream
&operator<<(std::ostream &os, const General_polygon_with_holes_2<Polygon_>& p)
{
  typename General_polygon_with_holes_2<Polygon_>::Hole_const_iterator hit;

  switch(get_mode(os)) {
    case IO::ASCII :
      os << p.outer_boundary() << ' ' << p.number_of_holes()<< ' ';
      for (hit = p.holes_begin(); hit != p.holes_end(); ++hit) {
        os << *hit << ' ';
      }
      return os;

    case IO::BINARY :
      os << p.outer_boundary()  << p.number_of_holes();
      for (hit = p.holes_begin(); hit != p.holes_end(); ++hit) {
        os << *hit;
      }
      return os;


    default:
      os << "General_polygon_with_holes_2( " << std::endl;
      os << p.outer_boundary() << " " << p.number_of_holes()<< " ";
      for (hit = p.holes_begin(); hit != p.holes_end(); ++hit) {
        os << *hit << " )";
      }
      return os;
  }
}

//-----------------------------------------------------------------------//
//                          operator>>
//-----------------------------------------------------------------------//

/*!
This operator imports a General_polygon_with_holes_2 from the input stream `in`.

An ASCII and a binary format exist. The stream detects the format
automatically and can read both.

The format consists of the number of curves of the outer boundary
followed by the curves themselves in counterclockwise order, followed
by the number of holes, and for each hole, the number of curves on its
outer boundary is followed by the curves themselves in clockwise
order.

\relates General_polygon_with_holes_2
*/
template <class Polygon_>
std::istream &operator>>(std::istream &is, General_polygon_with_holes_2<Polygon_>& p)
{
  p.clear();
  is >> p.outer_boundary();

  unsigned int n_holes;
  is >> n_holes;
  if (is)
  {
    Polygon_ pgn_hole;
    for (unsigned int i=0; i<n_holes; i++)
    {
      is >> pgn_hole;
      p.add_hole(pgn_hole);
    }
  }

  return is;
}


} //namespace CGAL

#endif
