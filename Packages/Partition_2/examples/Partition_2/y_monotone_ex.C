//------------------------------------------------------------------------------
//  y_monotone_ex
//
//    $Revision$
//    $Date$
//
//    program that computes a y-monotone partition of a particular polygon 
//    and checks that each polygon produced is, in fact, $y$-monotone and
//    that the polygons form a partition of the original polygon. 
//
//    (Note that the assertions are superfluous unless postcondition checking 
//    for the function y_monotone_partition_2 has been turned off.)
//------------------------------------------------------------------------------
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <list>

typedef CGAL::Cartesian<double>                           R;
typedef CGAL::Partition_traits_2<R>                       Traits;
typedef Traits::Point_2                                   Point_2;
typedef Traits::Polygon_2                                 Polygon_2;
typedef std::list<Polygon_2>                              Polygon_list;


void make_polygon(Polygon_2& polygon)
{
   polygon.push_back(Point_2(227,423));
   polygon.push_back(Point_2(123,364));
   polygon.push_back(Point_2(129,254));
   polygon.push_back(Point_2(230,285));
   polygon.push_back(Point_2(231,128));
   polygon.push_back(Point_2(387,205));
   polygon.push_back(Point_2(417,331));
   polygon.push_back(Point_2(319,225));
   polygon.push_back(Point_2(268,293));
   polygon.push_back(Point_2(367,399));
   polygon.push_back(Point_2(298,418));
   polygon.push_back(Point_2(196,326));
}


int main( )
{
   Polygon_2    polygon;
   Polygon_list partition_polys;

   make_polygon(polygon);
   CGAL::y_monotone_partition_2(polygon.vertices_begin(), 
                                polygon.vertices_end(),
                                std::back_inserter(partition_polys));

   std::list<Polygon_2>::const_iterator   poly_it;
   for (poly_it = partition_polys.begin(); poly_it != partition_polys.end();
        poly_it++)
   {
      assert(CGAL::is_y_monotone_2((*poly_it).vertices_begin(),
                                   (*poly_it).vertices_end()));
   }

   assert(CGAL::partition_is_valid_2(polygon.vertices_begin(),
                                     polygon.vertices_end(),
                                     partition_polys.begin(),
                                     partition_polys.end()));
   return 0;
}
