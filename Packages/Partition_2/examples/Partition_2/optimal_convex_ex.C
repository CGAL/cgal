//------------------------------------------------------------------------------
//  optimal_convex_ex
//
//    $Revision$
//    $Date$
//
//    program that computes an optimal convex partition of a particular 
//    polygon and checks the validity of the partition afterwards
//
//    (Note that the assertion is superfluous unless postcondition checking for
//    the function optimal_convex_partiton_2 has been turned off.)
//------------------------------------------------------------------------------
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/partition_2.h>
#include <list>

typedef CGAL::Cartesian<double>                           R;
typedef CGAL::Partition_traits_2<R>                       Traits;
typedef CGAL::Is_convex_2<Traits>                         Is_convex_2;
typedef Traits::Polygon_2                                 Polygon_2;
typedef Traits::Point_2                                   Point_2;
typedef Polygon_2::Vertex_const_iterator                  Vertex_iterator;
typedef std::list<Polygon_2>                              Polygon_list;
typedef CGAL::Partition_is_valid_traits_2<Traits, Is_convex_2>
                                                          Validity_traits;

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

int main()
{
   Polygon_2             polygon;
   Polygon_list          partition_polys;
   Traits                partition_traits;
   Validity_traits       validity_traits;

   make_polygon(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys),
                                    partition_traits);
   assert(CGAL::partition_is_valid_2(polygon.vertices_begin(), 
                                     polygon.vertices_end(),
                                     partition_polys.begin(), 
                                     partition_polys.end(),
                                     validity_traits));
   return 0;
}
