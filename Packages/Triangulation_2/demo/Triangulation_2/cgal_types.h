//CGAL headers
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>


typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>	    Rep;

typedef CGAL::Point_2<Rep>                  Point;
typedef CGAL::Segment_2<Rep>                Segment;
typedef CGAL::Line_2<Rep>                   Line;
typedef CGAL::Triangle_2<Rep>               Triangle;
typedef CGAL::Circle_2<Rep>                 Circle;

typedef CGAL::Triangulation_2<Rep>          Triangulation;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;


typedef Delaunay::Vertex_iterator           Vertex_iterator;
typedef Delaunay::Face_handle               Face_handle;
typedef Delaunay::Vertex_handle             Vertex_handle;
typedef Delaunay::Edge                      Edge;
typedef Delaunay::Line_face_circulator      Line_face_circulator;

