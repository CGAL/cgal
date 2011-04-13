//
//  file: examples/Partition_2/approx_convex_ex.C
//

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <cassert>
#include <list>

typedef double                                            NT;
typedef CGAL::Cartesian<NT>                               K;
typedef CGAL::Partition_traits_2<K>                       Traits;
typedef Traits::Point_2                                   Point_2;
typedef Traits::Polygon_2                                 Polygon_2;
typedef Polygon_2::Vertex_iterator                        Vertex_iterator;
typedef std::list<Polygon_2>                              Polygon_list;
typedef CGAL::Creator_uniform_2<int, Point_2>              Creator;
typedef CGAL::Random_points_in_square_2<Point_2, Creator> Point_generator;


int main()
{
   Polygon_2    polygon;
   Polygon_list partition_polys;

   CGAL::random_polygon_2(50, std::back_inserter(polygon), 
                          Point_generator(100));
   CGAL::approx_convex_partition_2(polygon.vertices_begin(), 
                                   polygon.vertices_end(),
                                   std::back_inserter(partition_polys));
   assert(CGAL::convex_partition_is_valid_2(polygon.vertices_begin(), 
                                            polygon.vertices_end(),
                                            partition_polys.begin(),
                                            partition_polys.end()));
   return 0;
}
