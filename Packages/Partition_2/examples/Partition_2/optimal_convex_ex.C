//
//  file: examples/Partition_2/optimal_convex_ex.C
//
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/partition_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <cassert>
#include <list>

typedef double                                            NT;
typedef CGAL::Cartesian<NT>                               K;
typedef CGAL::Partition_traits_2<K>                       Traits;
typedef CGAL::Is_convex_2<Traits>                         Is_convex_2;
typedef Traits::Polygon_2                                 Polygon_2;
typedef Traits::Point_2                                   Point_2;
typedef Polygon_2::Vertex_const_iterator                  Vertex_iterator;
typedef std::list<Polygon_2>                              Polygon_list;
typedef CGAL::Partition_is_valid_traits_2<Traits, Is_convex_2>
                                                          Validity_traits;
typedef CGAL::Creator_uniform_2<int, Point_2>              Creator;
typedef CGAL::Random_points_in_square_2<Point_2, Creator> Point_generator;

int main()
{
   Polygon_2             polygon;
   Polygon_list          partition_polys;
   Traits                partition_traits;
   Validity_traits       validity_traits;
   CGAL::random_polygon_2(50, std::back_inserter(polygon), 
                          Point_generator(100));
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
