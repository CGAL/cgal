#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/point_generators_2.h>


typedef CGAL::MP_Float				                      NT;
typedef CGAL::Cartesian<NT>                         K;
typedef CGAL::Partition_traits_2<K>                 Traits;
typedef Traits::Point_2                             Point_2;
typedef Traits::Polygon_2                           Cgal_Polygon;
typedef CGAL::Random_points_in_square_2<Point_2>    Point_generator;
