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

#ifndef CGAL_MIN_ELLIPSE_2_MIN_ELLIPSE_2_IMP_H
#define CGAL_MIN_ELLIPSE_2_MIN_ELLIPSE_2_IMP_H

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
              const Min_ellipse_2<Traits_>& min_ellipse)
{
    typedef typename Min_ellipse_2<Traits_>::Point  Point;
    typedef  std::ostream_iterator<Point>        Os_it;

    switch ( CGAL::IO::get_mode( os)) {

      case CGAL::IO::PRETTY:
        os << std::endl;
        os << "CGAL::Min_ellipse_2( |P| = "<<min_ellipse.number_of_points()
           << ", |S| = " << min_ellipse.number_of_support_points() << std::endl;
        os << "  P = {" << std::endl;
        os << "    ";
        std::copy( min_ellipse.points_begin(), min_ellipse.points_end(),
              Os_it( os, ",\n    "));
        os << "}" << std::endl;
        os << "  S = {" << std::endl;
        os << "    ";
        std::copy( min_ellipse.support_points_begin(),
              min_ellipse.support_points_end(),
              Os_it( os, ",\n    "));
        os << "}" << std::endl;
        os << "  ellipse = " << min_ellipse.ellipse() << std::endl;
        os << ")" << std::endl;
        break;

      case CGAL::IO::ASCII:
        std::copy( min_ellipse.points_begin(), min_ellipse.points_end(),
              Os_it( os, "\n"));
        break;

      case CGAL::IO::BINARY:
        std::copy( min_ellipse.points_begin(), min_ellipse.points_end(),
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
operator >> ( std::istream& is, CGAL::Min_ellipse_2<Traits_>& min_ellipse)
{
    switch ( CGAL::IO::get_mode( is)) {

      case CGAL::IO::PRETTY:
        std::cerr << std::endl;
        std::cerr << "Stream must be in ASCII or binary mode" << std::endl;
        break;

      case CGAL::IO::ASCII:
      case CGAL::IO::BINARY:
        typedef typename Min_ellipse_2<Traits_>::Point  Point;
        typedef  std::istream_iterator<Point>       Is_it;
        min_ellipse.clear();
        min_ellipse.insert( Is_it( is), Is_it());
        break;

      default:
        CGAL_assertion_msg( false, "CGAL::IO::mode invalid!");
        break; }

    return( is);
}

} //namespace CGAL

// ===== EOF ==================================================================

#endif // CGAL_MIN_ELLIPSE_2_MIN_ELLIPSE_2_IMP_H
