#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/point_generators_2.h>


typedef double                              Coord_type;
typedef CGAL::Cartesian<Coord_type>         Rep;

typedef CGAL::Point_2<Rep>                  Point;
typedef CGAL::Polygon_traits_2<Rep>         Traits;
typedef std::vector<Point>                  Container;
typedef CGAL::Polygon_2<Traits, Container>  Cgal_Polygon;
