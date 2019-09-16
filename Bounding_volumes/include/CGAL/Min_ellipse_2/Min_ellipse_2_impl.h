// Copyright (c) 1997-2001  
// ETH Zurich (Switzerland).  All rights reserved.
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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>, Bernd Gaertner

#include <iterator>

namespace CGAL {

// Class implementation (continued)
// ================================
// I/O
// ---
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os,
              const Min_ellipse_2<Traits_>& min_ellipse)
{
    using namespace std;

    typedef typename Min_ellipse_2<Traits_>::Point  Point;
    typedef  ostream_iterator<Point>        Os_it;

    switch ( CGAL::get_mode( os)) {

      case CGAL::IO::PRETTY:
        os << endl;
        os << "CGAL::Min_ellipse_2( |P| = "<<min_ellipse.number_of_points()
           << ", |S| = " << min_ellipse.number_of_support_points() << endl;
        os << "  P = {" << endl;
        os << "    ";
        copy( min_ellipse.points_begin(), min_ellipse.points_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  S = {" << endl;
        os << "    ";
        copy( min_ellipse.support_points_begin(),
              min_ellipse.support_points_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  ellipse = " << min_ellipse.ellipse() << endl;
        os << ")" << endl;
        break;

      case CGAL::IO::ASCII:
        copy( min_ellipse.points_begin(), min_ellipse.points_end(),
              Os_it( os, "\n"));
        break;

      case CGAL::IO::BINARY:
        copy( min_ellipse.points_begin(), min_ellipse.points_end(),
              Os_it( os));
        break;

      default:
        CGAL_optimisation_assertion_msg( false,
                                         "CGAL::get_mode( os) invalid!");
        break; }

    return( os);
}

template < class Traits_ >
std::istream&
operator >> ( std::istream& is, CGAL::Min_ellipse_2<Traits_>& min_ellipse)
{
    using namespace std;

    switch ( CGAL::get_mode( is)) {

      case CGAL::IO::PRETTY:
        cerr << endl;
        cerr << "Stream must be in ascii or binary mode" << endl;
        break;

      case CGAL::IO::ASCII:
      case CGAL::IO::BINARY:
        typedef typename Min_ellipse_2<Traits_>::Point  Point;
        typedef  istream_iterator<Point>       Is_it;
        min_ellipse.clear();
        min_ellipse.insert( Is_it( is), Is_it());
        break;

      default:
        CGAL_optimisation_assertion_msg( false, "CGAL::IO::mode invalid!");
        break; }

    return( is);
}

} //namespace CGAL

// ===== EOF ==================================================================
