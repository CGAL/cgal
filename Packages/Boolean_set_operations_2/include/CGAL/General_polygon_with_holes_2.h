// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef GENERAL_POLYGON_WITH_HOLES_2_H
#define GENERAL_POLYGON_WITH_HOLES_2_H

#include <list>

CGAL_BEGIN_NAMESPACE

template <class Polygon_>
class General_polygon_with_holes_2
{
public:

  typedef Polygon_                                     General_polygon_2;
    
  typedef std::list<General_polygon_2>                  Holes_containter;

  typedef typename Holes_containter::iterator          Holes_iterator;
  typedef typename Holes_containter::const_iterator    Holes_const_iterator;

  typedef unsigned int                                 Size;

  General_polygon_with_holes_2() : m_pgn()
  {}


  explicit General_polygon_with_holes_2(const General_polygon_2& pgn_boundary) : m_pgn(pgn_boundary)
  {}


  template <class HolesInputIterator>
  General_polygon_with_holes_2(const General_polygon_2& pgn_boundary,
                       HolesInputIterator h_begin,
                       HolesInputIterator h_end) : m_pgn(pgn_boundary),
                                                   m_holes(h_begin, h_end)
  {}

  Holes_iterator holes_begin()
  {
    return m_holes.begin();
  }

  Holes_iterator holes_end()
  {
    return m_holes.end();
  }

  Holes_const_iterator holes_begin() const
  {
    return m_holes.begin();
  }

  Holes_const_iterator holes_end() const
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

  void erase_hole(Holes_iterator hit)
  {
    m_holes.erase(hit);
  }

  bool has_holes() const
  {
    return (!m_holes.empty());
  }

  Size number_of_holes() const
  {
    return m_holes.size();
  }


protected:

  General_polygon_2           m_pgn;
  Holes_containter            m_holes;
};

CGAL_END_NAMESPACE

#endif
