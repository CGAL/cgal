// Copyright (c) 1997  ETH Zurich (Switzerland).
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
// Author(s)     : Sven Schoenherr <sven@inf.fu-berlin.de>
//                 Bernd Gaertner

#include <iterator>

namespace CGAL {

// Class implementation (continued)
// ================================
// I/O
// ---
template < class Traits >
std::ostream&
operator << ( std::ostream& os, const Min_sphere_d<Traits>& min_sphere)
{
    typedef typename Min_sphere_d<Traits>::Point  Point;

    switch ( get_mode( os)) {

      case IO::PRETTY:
        os << std::endl;
        os << "Min_sphere_d( |P| = " << min_sphere.number_of_points()
           << ", |S| = " << min_sphere.number_of_support_points()
           << std::endl;
        os << "  P = {" << std::endl;
        os << "    ";
        std::copy( min_sphere.points_begin(), min_sphere.points_end(),
              std::ostream_iterator<Point>( os, ",\n    "));
        os << "}" << std::endl;
        os << "  S = {" << std::endl;
        os << "    ";
        std::copy( min_sphere.support_points_begin(),
              min_sphere.support_points_end(),
              std::ostream_iterator<Point>( os, ",\n    "));
        os << "}" << std::endl;
        os << "  center = " << min_sphere.center() << std::endl;
        os << "  squared radius = " << min_sphere.squared_radius()
           << std::endl;
        os << ")" << std::endl;
        break;

      case IO::ASCII:
        os << min_sphere.number_of_points() << std::endl;
        std::copy( min_sphere.points_begin(), min_sphere.points_end(),
              std::ostream_iterator<Point>( os, "\n"));
        break;

      case IO::BINARY:
        os << min_sphere.number_of_points() << " ";
        std::copy( min_sphere.points_begin(), min_sphere.points_end(),
              std::ostream_iterator<Point>( os, " "));
        break;

      default:
        CGAL_optimisation_assertion_msg
            ( false, "get_mode( os) invalid!");
        break; }

    return( os);
}

template < class Traits >
std::istream&
operator >> ( std::istream& is, Min_sphere_d<Traits>& min_sphere)
{
    switch ( get_mode( is)) {

      case IO::PRETTY:
        std::cerr << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;

      case IO::ASCII:
      case IO::BINARY:
      {
        min_sphere.clear();
        int n; is >> n;
        typename Min_sphere_d<Traits>::Point p;
        for (int i=0; i<n; ++i) {
            is >> p;
            min_sphere.insert (p);
        }
      } break;

      default:
        CGAL_optimisation_assertion_msg( false, "IO::mode invalid!");
        break;
 }

    return( is);
}



} //namespace CGAL

// ===== EOF ==================================================================
