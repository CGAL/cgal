#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>
#include <list>

typedef CGAL::Cartesian<int>                        R;
typedef CGAL::Point_2<R>                            Point_2;
typedef CGAL::Polygon_2<R, std::list<Point_2> >     Polygon_2;

int main()
{
  Polygon_2                         CH;
  std::istream_iterator< Point_2 >  in_start( std::cin );
  std::istream_iterator< Point_2 >  in_end;

  CGAL::convex_hull_points_2( in_start, in_end,
                              std::inserter(CH, CH.vertices_begin()) );
  CGAL::Window_stream  W;
  W.init( -256.0, 255.0, -256.0);
  W.display();
  W << CH;

  W.read_mouse();          // wait for mouse click in window
  return 0;
}
