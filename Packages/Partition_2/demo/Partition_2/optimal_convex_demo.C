#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/IO/Window_stream.h>
#include <list>
#include <string>
#include <fstream>

typedef CGAL::Cartesian<double>                           R;
typedef CGAL::Partition_traits_2<R>                       Traits;
typedef Traits::Polygon_2                                 Polygon_2;
typedef Traits::Point_2                                   Point_2;
typedef Polygon_2::Vertex_const_iterator                  Vertex_iterator;
typedef std::list<Polygon_2>                              Polygon_list;

using CGAL::is_convex_2;

const int WINDOW_SIZE = 500;

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

void draw_polygons(const Polygon_2& polygon, 
                   const Polygon_list& partition_polys,
                   CGAL::Window_stream& W)
{
   W.clear();
   W.set_line_width(3);
   W << CGAL::BLACK;
   W << polygon; 

   W.set_line_width(1);
   W << CGAL::RED;
   Polygon_list::const_iterator p_it;
   for(p_it = partition_polys.begin();p_it != partition_polys.end(); p_it++)
      W << *p_it;
}


int main( int argc, char** argv )
{
   Polygon_2    polygon;
   Polygon_list partition_polys;
   Traits       partition_traits;
   CGAL::Window_stream W(WINDOW_SIZE, WINDOW_SIZE);

   W.init(0, WINDOW_SIZE, 0);

   CGAL::cgalize(W);

   W.display();

   make_polygon(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys),
                                    partition_traits);
   draw_polygons(polygon, partition_polys, W);
   // wait for mouse click
   W.read_mouse();
   return 0;
}
