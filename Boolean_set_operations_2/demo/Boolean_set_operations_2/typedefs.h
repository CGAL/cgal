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
#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

#include <CGAL/Cartesian.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/approximated_offset_2.h>

#ifdef CGAL_USE_GMP

  #include <CGAL/Gmpq.h>

  typedef CGAL::Gmpq                                    Base_nt;

#else

  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>

  typedef CGAL::Quotient<CGAL::MP_Float>                Base_nt;

#endif

typedef CGAL::Lazy_exact_nt<Base_nt>                  Coord_type;

// instead of
//typedef CGAL::Cartesian<Coord_type>		                Kernel;
// workaround for VC++
struct Kernel : public CGAL::Cartesian<Coord_type> {};

typedef Kernel::Segment_2			      Segment;
typedef Kernel::Point_2				      Point_2;
typedef Kernel::Circle_2                              Circle;
typedef Kernel::Iso_rectangle_2                       Iso_rectangle;

typedef CGAL::Gps_circle_segment_traits_2<Kernel>     Traits;
typedef Traits::Curve_2                               Curve;
typedef Traits::X_monotone_curve_2                    XCurve;
typedef Traits::Point_2                               Circular_point_2;
typedef Traits::Polygon_2                             Polygon_2;
typedef CGAL::General_polygon_with_holes_2<Polygon_2> Polygon_with_holes;
typedef CGAL::General_polygon_set_2<Traits>           Polygon_set;
typedef Polygon_with_holes::Hole_const_iterator       Hole_const_iterator;


typedef CGAL::Polygon_2<Kernel>                       Linear_polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>            Linear_polygon_with_holes_2;
#endif
