// Copyright (c) 1997-2001
// ETH Zurich (Switzerland).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>, Bernd Gaertner

#ifndef CGAL_MIN_CIRCLE_2_MIN_CIRCLE_2_IMPL_H
#define CGAL_MIN_CIRCLE_2_MIN_CIRCLE_2_IMPL_H

#include <CGAL/license/Bounding_volumes.h>

#include <iterator>

namespace CGAL {

// Class implementation (continued)
// ================================
// I/O
// ---
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os,
              const Min_circle_2<Traits_>& min_circle)
{
    typedef  typename Min_circle_2<Traits_>::Point  Point;
    typedef  std::ostream_iterator<Point>       Os_it;

    switch ( CGAL::IO::get_mode( os)) {

      case CGAL::IO::PRETTY:
        os << std::endl;
        os << "CGAL::Min_circle_2( |P| = " << min_circle.number_of_points()
           << ", |S| = " << min_circle.number_of_support_points() << std::endl;
        os << "  P = {" << std::endl;
        os << "    ";
        std::copy( min_circle.points_begin(), min_circle.points_end(),
              Os_it( os, ",\n    "));
        os << "}" << std::endl;
        os << "  S = {" << std::endl;
        os << "    ";
        std::copy( min_circle.support_points_begin(),
              min_circle.support_points_end(),
              Os_it( os, ",\n    "));
        os << "}" << std::endl;
        os << "  circle = " << min_circle.circle() << std::endl;
        os << ")" << std::endl;
        break;

      case CGAL::IO::ASCII:
        std::copy( min_circle.points_begin(), min_circle.points_end(),
              Os_it( os, "\n"));
        break;

      case CGAL::IO::BINARY:
        std::copy( min_circle.points_begin(), min_circle.points_end(),
              Os_it( os));
        break;

      default:
        CGAL_assertion_msg( false,
                                         "CGAL::IO::get_mode( os) invalid!");
        break; }

    return( os);
}

template < class Traits_ >
std::istream&
operator >> ( std::istream& is, CGAL::Min_circle_2<Traits_>& min_circle)
{
    switch ( CGAL::IO::get_mode( is)) {

      case CGAL::IO::PRETTY:
        std::cerr << std::endl;
        std::cerr << "Stream must be in ASCII or binary mode" << std::endl;
        break;

      case CGAL::IO::ASCII:
      case CGAL::IO::BINARY:
        typedef  typename CGAL::Min_circle_2<Traits_>::Point  Point;
        typedef  std::istream_iterator<Point>            Is_it;
        min_circle.clear();
        min_circle.insert( Is_it( is), Is_it());
        break;

      default:
        CGAL_assertion_msg( false, "CGAL::IO::mode invalid!");
        break; }

    return( is);
}

} //namespace CGAL

// ===== EOF ==================================================================

#endif // CGAL_MIN_CIRCLE_2_MIN_CIRCLE_2_IMPL_H
