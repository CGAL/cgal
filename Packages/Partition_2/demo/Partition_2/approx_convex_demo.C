
// file: demo/Partition_2/approx_convex_demo.C
//
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/IO/Window_stream.h>
#include <list>
#include <string>
#include <fstream>

typedef CGAL::Cartesian<double>                           K;
typedef CGAL::Polygon_traits_2<K>                         Traits;
typedef Traits::Point_2                                   Point_2;
typedef std::list<Point_2>                                Container;
typedef CGAL::Polygon_2<Traits, Container>                Polygon_2;
typedef Polygon_2::Vertex_iterator                        Vertex_iterator;
typedef std::list<Polygon_2>                              Polygon_list;

const int WINDOW_SIZE = 500;
enum Button_nums {EXIT = 4, REFRESH};  // button numbers for window's buttons

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
   for(p_it = partition_polys.begin(); p_it != partition_polys.end(); p_it++)
      W << *p_it;
}


int main( )
{
   Polygon_2    polygon;
   Polygon_list partition_polys;
   CGAL::Window_stream W(WINDOW_SIZE, WINDOW_SIZE, 
                         "Approximately Optimal Convex Partition");

   W.init(0, WINDOW_SIZE, 0);
   W.button("Refresh", REFRESH);
   W.button("Exit", EXIT);

   W.display();

   make_polygon(polygon);
   CGAL::approx_convex_partition_2(polygon.vertices_begin(), 
                                   polygon.vertices_end(),
                                   std::back_inserter(partition_polys));
   while (1) // wait for click on exit button
   {
      draw_polygons(polygon, partition_polys, W);
      if (W.read_mouse() == EXIT)
        return 0;
   }
}

