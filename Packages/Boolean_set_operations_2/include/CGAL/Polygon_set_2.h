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

#ifndef POLYGON_SET_2_H
#define POLYGON_SET_2_H

CGAL_BEGIN_NAMESPACE

#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <verctor>

template <class Kernel,
          class Containter = std:vector<typename Kernel::Point_2> >
class Polygon_set_2 :
  public General_polygon_set_2<Gps_segment_traits_2<Kernel, Containter> >
{
  typedef General_polygon_set_2<Gps_segment_traits_2<Kernel, Containter> > Base;
  typedef Bso_dcel<Traits_2>                                Bso_dcel;  
  
public:

  typedef Arrangement_2<Traits_2, Bso_dcel>                     Arrangement_2;
  typedef Gps_segment_traits_2<Kernel, Containter>::Traits_2    Traits_2;
  typedef typename Traits_2::Polygon_2                          Polygon_2;
  typedef typename Traits_2::Polygon_with_holes_2               Polygon_with_holes_2;

  Polygon_set_2() : Base()
  {}

  // constructor with traits object
  Polygon_set_2(Traits_2& tr) : Base(tr)
  {}

  explicit Polygon_set_2(const Polygon_2& pgn) : Base(pgn)
  {}

  explicit Polygon_set_2(const Polygon_with_holes_2& pgn): Base(pgn)
  {}

  template <class PolygonIterator>
  Polygon_set_2(PolygonIterator pgn_begin, PolygonIterator pgn_end):
    Base(pgn_begin, pgn_end)
  {} 

  template <class PolygonIterator, class PolygonWithHolesIterator>
  Polygon_set_2(PolygonIterator pgn_begin,
                PolygonIterator pgn_end,
                PolygonWithHolesIterator  pgn_with_holes_begin,
                PolygonWithHolesIterator  pgn_with_holes_end):
    Base(pgn_begin, pgn_end, pgn_with_holes_begin, pgn_with_holes_end)
  {}

};

CGAL_END_NAMESPACE

#endif
