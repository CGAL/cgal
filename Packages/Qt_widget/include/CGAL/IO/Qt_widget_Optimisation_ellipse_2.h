// Copyright (c) 1997-2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
        ws << oe.to_double();
        break;
      default:
        CGAL_optimisation_assertion( ( oe.n_boundary_points >= 0) &&
                                     ( oe.n_boundary_points <= 5) ); }
    return( ws);
}

}//end namespace CGAL

#endif





