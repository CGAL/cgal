// tell the traits class to provide special functionality ...
#define CGAL_PROVIDE_LEDA_PARTITION_TRAITS
// turn the check off !!!
#define CGAL_NO_POSTCONDITIONS


#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CEP/Leda_rat_kernel/polygon_Window_stream_spec.h>
#include <CGAL/partition_2.h>

#include <LEDA/rat_window.h> 
#include <CGAL/IO/Window_stream.h>

#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <list>
#include <string>
#include <fstream>

typedef CGAL::leda_rat_kernel_traits                      K;
typedef K::Point_2                                        Point_2;
typedef K::Polygon_2                                      Polygon_2;
typedef Polygon_2::Vertex_iterator                        Vertex_iterator;
typedef std::list<Polygon_2>                              Polygon_list;

const int WINDOW_SIZE = 500;
enum Button_nums {EXIT = 4, REFRESH};  // button numbers for window's buttons

void make_polygon(Polygon_2& polygon)
{
   polygon.push_back(Point_2(227,423,1));
   polygon.push_back(Point_2(123,364,1));
   polygon.push_back(Point_2(129,254,1));
   polygon.push_back(Point_2(230,285,1));
   polygon.push_back(Point_2(231,128,1));
   polygon.push_back(Point_2(387,205,1));
   polygon.push_back(Point_2(417,331,1));
   polygon.push_back(Point_2(319,225,1));
   polygon.push_back(Point_2(268,293,1));
   polygon.push_back(Point_2(367,399,1));
   polygon.push_back(Point_2(298,418,1));
   polygon.push_back(Point_2(196,326,1));
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


int main( )
{
   Polygon_2    polygon;
   Polygon_list partition_polys;
   CGAL::Window_stream W(WINDOW_SIZE, WINDOW_SIZE,
                         "Y-Monotone Partition");

   W.init(0, WINDOW_SIZE, 0);
   W.button("Refresh", REFRESH);
   W.button("Exit", EXIT);

   W.display();

   make_polygon(polygon);
   
   K leda_traits;
   
   CGAL::approx_convex_partition_2(polygon.vertices_begin(), 
                                polygon.vertices_end(),
                                std::back_inserter(partition_polys),
				leda_traits);
   draw_polygons(polygon, partition_polys, W);

   while (1) // wait for click on exit button
   {
      draw_polygons(polygon, partition_polys, W);
      if (W.read_mouse() == EXIT)
        return 0;
   }

   return 0;
}

