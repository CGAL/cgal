#include <iostream>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

template<class K>
void read_polygon( std::istream& is, CGAL::Orientation aO, CGAL::Polygon_2<K>& rPoly )
{
  is >> rPoly ;
  if ( rPoly.orientation() != aO )
  {
    rPoly.reverse_orientation();
  }
}

template<class K>
void read_polygon_with_holes( std::istream& is, CGAL::Polygon_with_holes_2<K>& rPolyWH )
{
   int ccb_count = 0 ;
   is >> ccb_count ;
   read_polygon( is, CGAL::COUNTERCLOCKWISE, rPolyWH.outer_boundary() ) ;
   for ( int i = 1 ; i < ccb_count ; ++ i )
   {
     CGAL::Polygon_2<K> lHole ;
     read_polygon( is, CGAL::CLOCKWISE, lHole ) ;
     rPolyWH.add_hole(lHole);
   }
}

