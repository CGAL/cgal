#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/point_generators_2.h>


typedef double                                   Coord_type;
typedef CGAL::Cartesian<Coord_type>              K;
typedef K::Point_2                               Point_2;
typedef std::vector<Point_2>                     Container;
typedef CGAL::Polygon_2<K, Container>            Cgal_Polygon;
typedef CGAL::Random_points_in_square_2<Point_2> Point_generator;
