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

#ifndef CGAL_POLYGON_WITH_HOLES_2_H
#define CGAL_POLYGON_WITH_HOLES_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class Kernel,
          class Containter = std::vector<typename Kernel::Point_2> >
class Polygon_with_holes_2 : 
  public General_polygon_with_holes_2<CGAL::Polygon_2<Kernel, Containter> >
{
public:

  typedef CGAL::Polygon_2<Kernel, Containter>        Polygon_2;
  typedef General_polygon_with_holes_2<Polygon_2>    Base;
  typedef typename Base::Hole_const_iterator         Hole_const_iterator;
  typedef typename Base::Size                        Size;

  /*! Default constructor. */
  Polygon_with_holes_2 () : 
    Base()
  {}

  /*! Constructor from the base class. */
  Polygon_with_holes_2 (const Base& base) : 
    Base (base)
  {}

  /*! Constructor from a polygon. */
  explicit Polygon_with_holes_2 (const Polygon_2& pgn_boundary) : 
    Base (pgn_boundary)
  {}

  /*! Constructor from a polygon (outer boundary) and hole polygons. */
  template <class HolesInputIterator>
  Polygon_with_holes_2 (const Polygon_2& pgn_boundary,
                        HolesInputIterator h_begin,
                        HolesInputIterator h_end) : 
    Base (pgn_boundary, h_begin, h_end)
  {}

};

CGAL_END_NAMESPACE

#endif
