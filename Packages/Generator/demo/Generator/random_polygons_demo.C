// file : demo/Generator/random_polygons_demo.C
// --------------------------------------------
// program generating random simple polygons

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>
#include <fstream>

#if !defined(CGAL_USE_CGAL_WINDOW)
#include <LEDA/menu.h>
#include <LEDA/panel.h>
#include <LEDA/point.h>
#include <LEDA/string.h>
#else
#include <CGAL/LEDA/menu.h>
#include <CGAL/LEDA/panel.h>
#include <string>
#endif

#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/copy_n.h>
#include <CGAL/Aff_transformation_2.h>
#include <list>
#include <vector>

#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_panel CGAL::panel
#endif

typedef double                             NT;
typedef CGAL::Cartesian<NT>                K;
typedef K::Vector_2                        Vector_2;
typedef CGAL::Polygon_traits_2<K>          Traits;
typedef Traits::Point_2                    Point_2;
typedef std::list<Point_2>                 Container;
typedef CGAL::Polygon_2<Traits, Container> Polygon_2;
typedef Polygon_2::Vertex_iterator         Vertex_iterator;
typedef Polygon_2::Vertex_const_iterator   Vertex_const_iterator;
typedef CGAL::Aff_transformation_2<K>      Transformation_2;

enum Button_nums {QUIT=4, RANDOM_POINT_SET, INPUT_POINT_SET, SIMPLE_POLYGON};

const int WINDOW_SIZE = 500;
const int MAX_POLY_SIZE = 100;
const double RADIUS = 0.9 * WINDOW_SIZE/2.0;

Transformation_2 translate(CGAL::TRANSLATION, 
                           Vector_2(WINDOW_SIZE/2.0, WINDOW_SIZE/2.0));

typedef CGAL::Creator_uniform_2<double, Point_2>             Creator;
typedef CGAL::Random_points_in_square_2<Point_2, Creator> Point_generator;

CGAL::Window_stream W(WINDOW_SIZE, WINDOW_SIZE, "Random polygons");
leda_panel controls("Controls");

template <class OutputIterator>
OutputIterator input_point_set(OutputIterator result)
{
   W.acknowledge("Click right mouse button to terminate input of points");
   W << CGAL::RED;
   double x, y;
   while (W.read_mouse(x, y) != MOUSE_BUTTON(3))
   {
      Point_2 point(x,y);
      W << point;
      *result = point;
      result++;
   }
   return result;
}


template <class ForwardIterator>
void draw_points(ForwardIterator first, ForwardIterator beyond, 
                 CGAL::Color color)
{
   int i = 0;
   W << color;
   for (ForwardIterator it = first; it != beyond; it++, i++)
   {
#if !defined(CGAL_USE_CGAL_WINDOW)   
       leda_point point(CGAL::to_double((*it).x()),CGAL::to_double((*it).y()));
       leda_string label("%d", i);
       W.draw_text(point, label);
#endif
       W << *it;
   }
}

void draw_polygon(const Polygon_2& polygon, CGAL::Color color, int line_width)
{
   W.set_line_width(line_width);
   W << color;
   W << polygon;
   draw_points(polygon.vertices_begin(), polygon.vertices_end(), CGAL::RED);
}


int main( )
{
   typedef std::vector<Point_2>   Point_vector;
   typedef Point_vector::iterator Point_iterator;

   Polygon_2             polygon;
   Point_vector          point_set;

   W.init(0, WINDOW_SIZE, 0);
   W.set_node_width(3);
   int size = 4;

   controls.int_item("Number of vertices", size, 4, MAX_POLY_SIZE);
   controls.button("Random Points", RANDOM_POINT_SET);
   controls.button("Simplify", SIMPLE_POLYGON);
   controls.button("Input Points", INPUT_POINT_SET);
   controls.button("Quit", QUIT);

   int button_num;
   bool polygon_computed = false;

   W.display();
   controls.display();
   while (1)
   {
      W.clear();
      if (polygon_computed)
         draw_polygon(polygon, CGAL::BLACK, 2);
      else
         draw_points(point_set.begin(), point_set.end(), CGAL::RED);
      button_num = controls.read_mouse();
      switch (button_num) 
      {
         case QUIT:
            return(0);
         case INPUT_POINT_SET:
         {
            point_set.clear();
            input_point_set(std::back_inserter(point_set));
            polygon_computed = false;
            break;
         }
         case RANDOM_POINT_SET:
         {
            point_set.clear();
            CGAL::copy_n(Point_generator(RADIUS), size,
                         std::back_inserter(point_set)); 
            Point_iterator it;
            for (it = point_set.begin(); it != point_set.end(); it++)
               *it = translate.transform(*it);
            polygon_computed = false;
            break;
         }
         case SIMPLE_POLYGON:
         {
            if (point_set.empty()) break;

            if (!CGAL::duplicate_points(point_set.begin(), point_set.end()) )
            {
               std::random_shuffle(point_set.begin(), point_set.end());
               polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
               CGAL::random_polygon_2(point_set.size(), 
                                      std::back_inserter(polygon), 
                                      point_set.begin());
               polygon_computed = true;
            }
            else
            {
               W.acknowledge("Duplicate points; simplification not possible");
            }
            break;
         }
         default:
            break;
      }
    }
    return 0; // statement unreachable; warning ok
}

