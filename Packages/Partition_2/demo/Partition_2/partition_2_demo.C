//
// file: demo/Partition_2/partition_2_demo.C
//
#include <CGAL/basic.h>
#ifndef CGAL_USE_LEDA
int main() { 
   std::cout << "\nSorry, this demo needs LEDA for visualization.\n"; 
   return 0; 
}
#else
#include <typedefs.h>
#include <polygon_io.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/copy_n.h>
#include <CGAL/Random.h>
#include <CGAL/random_polygon_2.h>
#include <LEDA/menu.h>
#include <LEDA/string.h>
#include <LEDA/polygon.h>
#include <LEDA/point.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>
#include <fstream>

enum Button_nums {QUIT=4, MOUSE_POLYGON, RANDOM_POLYGON, 
                  READ_FROM_FILE, WRITE_INPUT_TO_FILE};

const int WINDOW_SIZE = 500;
const int MAX_POLY_SIZE = 100;
const double RADIUS = 0.9 * WINDOW_SIZE/2.0;

Transformation_2 translate(CGAL::TRANSLATION,
                           Vector_2(WINDOW_SIZE/2.0, WINDOW_SIZE/2.0));

Polygon_2 polygon;
std::list<Polygon_2> greene_approx_polys;
std::list<Polygon_2> greene_opt_polys;
std::list<Polygon_2> tri_approx_polys;
std::list<Polygon_2> monotone_polys;
CGAL::Window_stream W(WINDOW_SIZE, WINDOW_SIZE, "Polygon decomposition");
panel_item file_name_panel;
leda_string file_label("Current file: ");
leda_string current_file_name("data/_last");
leda_panel controls(WINDOW_SIZE, 10, "Controls");

bool show_greene_approx = false;
bool show_monotone      = false;
bool show_greene_opt    = false;
bool show_tri_approx    = false;
bool show_coords        = false;


// button_num is necessary for using this as an action function for the
// bool_item panel buttons
void draw_polygons(int button_num)
{
   Polygon_2 poly;
   unsigned int width = static_cast<unsigned int>(show_greene_approx) +
                        static_cast<unsigned int>(show_monotone) +
                        static_cast<unsigned int>(show_greene_opt) +
                        static_cast<unsigned int>(show_tri_approx);

   double legend_dist = 25;
   double legend_x = W.xmin(); 
   double legend_y = W.ymin();

   W.clear();
   W.set_line_width(width+2);
   W << CGAL::BLACK;
   if (!polygon.is_empty())  
   {
      draw_a_polygon(W, polygon, CGAL::BLACK);
      width--;
      if (show_coords) label_vertices(W, polygon);
   }
   legend_y += legend_dist;

   int num_greene_opt;
   if (show_greene_opt) 
   {
      draw_poly_list(W, greene_opt_polys, width, CGAL::GREEN, num_greene_opt,
                     show_coords);
      leda_string go_num_label(" greene opt. (%d)", num_greene_opt);
      W.draw_text(legend_x, legend_y, go_num_label);
      legend_y += legend_dist;
      width--;
   }
   int num_tri_approx;
   if (show_tri_approx) 
   {
      draw_poly_list(W, tri_approx_polys, width, CGAL::VIOLET, num_tri_approx,
                     show_coords);
      leda_string tri_num_label(" tri. approx. (%d)", num_tri_approx);
      W.draw_text(legend_x, legend_y, tri_num_label);
      legend_y += legend_dist;
      width--;
   }
   int num_greene_approx;
   if (show_greene_approx) 
   {
      draw_poly_list(W, greene_approx_polys, width, CGAL::BLUE, 
                     num_greene_approx, show_coords);
      leda_string ga_num_label(" greene approx. (%d)", num_greene_approx);
      W.draw_text(legend_x, legend_y, ga_num_label);
      legend_y += legend_dist;
      width--;
   }
   int num_monotone; 
   if (show_monotone) 
   {
      draw_poly_list(W, monotone_polys, width, CGAL::RED, num_monotone,
                     show_coords);
      leda_string mono_num_label(" monotone (%d)", num_monotone);
      W.draw_text(legend_x, legend_y, mono_num_label);
      legend_y += legend_dist;
      width--;
   }

}


void compute_partitions(Polygon_2& polygon)
{
    CGAL::Partition_traits_2<K> partition_traits;
    save_poly_to_file(polygon, "data/_last");
    monotone_polys.clear();
    greene_approx_polys.clear();
    greene_opt_polys.clear();
    tri_approx_polys.clear();

    if (polygon.is_empty()) return;

    if (CGAL::orientation_2(polygon.vertices_begin(), polygon.vertices_end(), 
                      partition_traits) != CGAL::COUNTERCLOCKWISE)
       std::reverse(polygon.vertices_begin(), polygon.vertices_end());


    CGAL::approx_convex_partition_2(polygon.vertices_begin(), 
                                    polygon.vertices_end(), 
                                    std::back_inserter(tri_approx_polys), 
                                    partition_traits);
    save_partition_p_list_to_file(tri_approx_polys, current_file_name,
                                  ".approx_cvx");

    CGAL::y_monotone_partition_2(polygon.vertices_begin(), 
                                 polygon.vertices_end(),
                                 std::back_inserter(monotone_polys), 
                                 partition_traits);
    save_partition_p_list_to_file(monotone_polys, current_file_name,
                                  ".y_mono");

    CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                     polygon.vertices_end(),
                                     std::back_inserter(greene_opt_polys), 
                                     partition_traits);
    save_partition_p_list_to_file(greene_opt_polys, current_file_name,
                                  ".opt_cvx");

    CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(), 
                                        polygon.vertices_end(),
                                        std::back_inserter(greene_approx_polys),
                                        partition_traits);
    save_partition_p_list_to_file(greene_opt_polys, current_file_name,
                                  ".greene_approx");

}

int main( )
{
   menu polygon_menu;
   polygon_menu.button("Enter new polygon", MOUSE_POLYGON);
   polygon_menu.button("Random polygon", RANDOM_POLYGON);

   W.init(0, WINDOW_SIZE, 0);

   controls.button("Load", READ_FROM_FILE);
   controls.button("Save Input", WRITE_INPUT_TO_FILE);
   controls.button("Polygon", polygon_menu);
   controls.button("Quit", QUIT);
   controls.make_menu_bar();
   controls.bool_item("Show convex optimal", show_greene_opt, draw_polygons);
   controls.bool_item("Show convex approximation", show_tri_approx, 
                      draw_polygons);
   controls.bool_item("Show Greene approximation", show_greene_approx, 
                      draw_polygons);
   controls.bool_item("Show monotone partition", show_monotone, draw_polygons);
   controls.bool_item("Show coordinates", show_coords, draw_polygons);
   file_name_panel =
       controls.text_item(file_label.insert(file_label.length(), 
                          current_file_name));

   W.display();
   controls.display();
   int button_num;

   if (read_poly_from_file(polygon, current_file_name))
      compute_partitions(polygon);

   W.display();
   while (1)
   {
      draw_polygons(button_num);
      button_num = controls.read_mouse();
      switch (button_num) 
      {
         case QUIT:
            return(0);
         case MOUSE_POLYGON:
         {
            input_polygon(polygon, W);
            compute_partitions(polygon);
            break;
         }
         case RANDOM_POLYGON:
         {
            polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
            CGAL::random_polygon_2(CGAL::Random().get_int(4,MAX_POLY_SIZE),
                                   std::back_inserter(polygon),
                                   Point_generator(RADIUS));
            polygon = CGAL::transform(translate, polygon);
            compute_partitions(polygon);
            break;
         }
         case READ_FROM_FILE: 
         {
            leda_string new_file_name = W.read_string("Load from file:");
            if (read_poly_from_file(polygon, new_file_name)) 
            {
               compute_partitions(polygon);
               current_file_name = new_file_name;
               controls.set_text(file_name_panel, 
                              file_label.insert(file_label.length(),
                                                current_file_name));
               controls.redraw_panel();
            }
            else
            {
               W.acknowledge(new_file_name + ": No such file or directory");
            }
            break;
         }
         case WRITE_INPUT_TO_FILE:
         {
            current_file_name = W.read_string("Save input to file:");
            controls.set_text(file_name_panel, 
                              file_label.insert(file_label.length(),
                                                current_file_name));
            controls.redraw_panel();
            save_poly_to_file(polygon, current_file_name);
            break;
         }
         default:
            break;
      }
    }
    return 0; // statement is unreachable; warning is OK
}

#endif // CGAL_USE_LEDA
