// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Radu Ursu

#ifndef CGAL_QT_WIDGET_OPTIMISATION_ELLIPSE_2_H
#define CGAL_QT_WIDGET_OPTIMISATION_ELLIPSE_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Conic_2.h>

namespace CGAL{

template< class Traits_ >
Qt_widget&
operator << ( Qt_widget &ws,
              const CGAL::Optimisation_ellipse_2<Traits_>& oe)
{

  typedef Cartesian<double> Rep;
  typedef Point_2<Rep>	    Point;
  typedef Segment_2<Rep>    Segment;
  
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
        ws << Segment( Point(px1, py1), Point(px2, py2)); 
	      }
        break;
      case 3:
      case 4:
      case 5:
	{
	  typedef CGAL::Conic_2<CGAL::Cartesian<double> >
	    DoubleConic_2;
	  DoubleConic_2 dc2;
	  oe.double_conic(dc2);
	  ws << dc2;
	}
        break;
      default:
        CGAL_optimisation_assertion( ( oe.n_boundary_points >= 0) &&
                                     ( oe.n_boundary_points <= 5) ); }
    return( ws);
}

}//end namespace CGAL

#endif
