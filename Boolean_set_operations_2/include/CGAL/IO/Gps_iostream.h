// Copyright (c) 2009  Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_GPS_IOSTREAM_H
#define CGAL_GPS_IOSTREAM_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

#include <iostream>
#include <list>

#include <CGAL/basic.h>
#include <CGAL/General_polygon_set_2.h>

namespace CGAL {

template <typename Traits>
std::ostream & operator<< (std::ostream& os,
                           const CGAL::General_polygon_set_2<Traits> & pgn_set)
{
  typedef typename CGAL::General_polygon_set_2<Traits>::Polygon_with_holes_2
                                                        Polygon_with_holes_2;
  typedef std::list<Polygon_with_holes_2>               Pgn_with_holes_container;

  Pgn_with_holes_container res;
  pgn_set.polygons_with_holes (std::back_inserter (res));

  std::cout << pgn_set.number_of_polygons_with_holes() << std::endl;
  std::copy(res.begin(), res.end(),
            std::ostream_iterator<Polygon_with_holes_2>(std::cout, "\n"));
  
  return os;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
