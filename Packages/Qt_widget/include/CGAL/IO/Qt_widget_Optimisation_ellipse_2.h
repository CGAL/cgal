#ifndef CGAL_QT_WIDGET_OPTIMISATION_ELLIPSE_2_H
#define CGAL_QT_WIDGET_OPTIMISATION_ELLIPSE_2_H

#include <CGAL/IO/Qt_widget>
#include <CGAL/Optimisation_ellipse_2.h>

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





