#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Largest_empty_iso_rectangle_2.h>

typedef double                      Coord_type;
typedef CGAL::Cartesian<Coord_type> Rep;

typedef Rep::Point_2             Point;
typedef Rep::Segment_2           Segment;
typedef Rep::Iso_rectangle_2     Iso_rectangle_2;
typedef CGAL::Largest_empty_iso_rectangle_2<Rep>
                                 Largest_empty_rect;
