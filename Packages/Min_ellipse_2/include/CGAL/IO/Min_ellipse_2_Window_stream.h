// Copyright (c) 1997-2001  Freie Universitaet Berlin (Germany).
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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>, Bernd Gaertner

// Each of the following operators is individually
// protected against multiple inclusion.

// Window_stream I/O operators
// ===========================
// includes
#include <CGAL/IO/Conic_2_Window_stream.h>
#include <CGAL/Cartesian.h>
// Optimisation_ellipse_2
// ----------------------
#ifdef CGAL_OPTIMISATION_ELLIPSE_2_H
#ifndef CGAL_IO_WINDOW_STREAM_OPTIMISATION_ELLIPSE_2
#define CGAL_IO_WINDOW_STREAM_OPTIMISATION_ELLIPSE_2

template< class Traits_ >
CGAL::Window_stream&
operator << ( CGAL::Window_stream &ws,
              const CGAL::Optimisation_ellipse_2<Traits_>& oe)
{
    switch ( oe.n_boundary_points) {
      case 0:
        break;
      case 1:
        ws << oe.boundary_point1;
        break;
      case 2: {
        double  px1( CGAL::to_double( oe.boundary_point1.x()));
        double  py1( CGAL::to_double( oe.boundary_point1.y()));
        double  px2( CGAL::to_double( oe.boundary_point2.x()));
        double  py2( CGAL::to_double( oe.boundary_point2.y()));
        ws.draw_segment( px1, py1, px2, py2); }
        break;
      case 3:
      case 4:
      case 5:
	Cartesian<double>::Conic_2 dc;
	oe.double_conic(dc);
        ws << dc;
        break;
      default:
        CGAL_optimisation_assertion( ( oe.n_boundary_points >= 0) &&
                                     ( oe.n_boundary_points <= 5) ); }
    return( ws);
}

#endif // CGAL_IO_WINDOW_STREAM_OPTIMISATION_ELLIPSE_2
#endif // CGAL_OPTIMISATION_ELLIPSE_2_H

// Min_ellipse_2
// -------------
#ifdef CGAL_MIN_ELLIPSE_2_H
#ifndef CGAL_IO_WINDOW_STREAM_MIN_ELLIPSE_2
#define CGAL_IO_WINDOW_STREAM_MIN_ELLIPSE_2

template< class Traits_ >
CGAL::Window_stream&
operator << ( CGAL::Window_stream &ws,
              const CGAL::Min_ellipse_2<Traits_>& min_ellipse)
{
    typedef typename CGAL::Min_ellipse_2<Traits_>::Point_iterator
	    Point_iterator;

    Point_iterator  first( min_ellipse.points_begin());
    Point_iterator  last ( min_ellipse.points_end());
    for ( ; first != last; ++first)
        ws << *first;
    return( ws << min_ellipse.ellipse());
}
#endif // CGAL_IO_WINDOW_STREAM_MIN_ELLIPSE_2
#endif // CGAL_MIN_ELLIPSE_2_H

// ===== EOF ==================================================================
