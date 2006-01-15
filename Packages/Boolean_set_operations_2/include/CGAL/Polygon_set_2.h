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

#ifndef CGAL_POLYGON_SET_2_H
#define CGAL_POLYGON_SET_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class Kernel,
          class Containter = std::vector<typename Kernel::Point_2> >
class Polygon_set_2 :
  public General_polygon_set_2<Gps_segment_traits_2<Kernel, Containter> >
{  
private:
  typedef General_polygon_set_2<Gps_segment_traits_2<Kernel, Containter> >                 
                                                          Base;

public:
  typedef  typename Base::Traits_2                        Traits_2;    
  typedef  typename Base::Polygon_2                       Polygon_2;
  typedef  typename Base::Polygon_with_holes_2            Polygon_with_holes_2;
  typedef  typename Base::Arrangement_2                   Arrangement_2;
  typedef  typename Base::Size                            Size;

  /*! Default consturctor. */
  Polygon_set_2 () :
    Base()
  {}

  /*! Consturctor from the base class. */
  Polygon_set_2 (const Base& base) :
    Base (base)
  {}

  /*! Constructor with traits object. */
  Polygon_set_2 (Traits_2& tr) :
    Base(tr)
  {}

  /*! Constructor from a polygon. */
  explicit Polygon_set_2 (const Polygon_2& pgn) :
    Base (pgn)
  {}

  /*! Constructor from a polygon with holes. */
  explicit Polygon_set_2 (const Polygon_with_holes_2& pwh):
    Base (pwh)
  {}

  /*! Constructor from a range of polygons. */
  template <class PolygonIterator>
  Polygon_set_2 (PolygonIterator pgn_begin, PolygonIterator pgn_end) :
    Base(pgn_begin, pgn_end)
  {} 

  /*! Constructor from ranges of polygons and polygons with holes. */
  template <class PolygonIterator, class PolygonWithHolesIterator>
  Polygon_set_2 (PolygonIterator pgn_begin,
                 PolygonIterator pgn_end,
                 PolygonWithHolesIterator pgn_with_holes_begin,
                 PolygonWithHolesIterator pgn_with_holes_end):
    Base (pgn_begin, pgn_end, 
          pgn_with_holes_begin, pgn_with_holes_end)
  {}

};

CGAL_END_NAMESPACE

#endif
