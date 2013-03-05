// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

namespace CGAL {

template <class Polygon_>
class General_polygon_with_holes_2
{

public:
  typedef General_polygon_with_holes_2<Polygon_> 		Self;
  typedef Polygon_                                    Polygon_2;
  typedef typename Self::Polygon_2							General_polygon_2;  
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

template <class Polygon_>
std::ostream
&operator<<(std::ostream &os, const General_polygon_with_holes_2<Polygon_>& p)
{
  typename General_polygon_with_holes_2<Polygon_>::Hole_const_iterator hit;

  switch(os.iword(IO::mode)) {
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
