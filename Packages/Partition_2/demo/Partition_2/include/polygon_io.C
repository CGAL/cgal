#include <polygon_io.h>


void convert_leda_polygon_to_CGAL(const leda_polygon &ledaP, 
                                  Polygon_2 &p)
{
   leda_point v;
   p.erase(p.vertices_begin(), p.vertices_end());
   forall_vertices(v,ledaP)
   {
      p.push_back(Point_2(v.xcoord(),v.ycoord()));
   }
}

void convert_point_list_to_leda_poly(const Polygon_2& p,leda_polygon &ledaP) 
{
   Polygon_2::const_iterator it;
   leda_list<leda_point> L;
   for(it = p.vertices_begin(); it != p.vertices_end(); it++)
   {
      L.push_back(leda_point(CGAL::to_double((*it).x()),
                             CGAL::to_double((*it).y())));
   }
   leda_polygon temp(L);
   ledaP = temp;
}

void input_polygon(Polygon_2& p, CGAL::Window_stream& W)
{
    W << CGAL::BLACK;
    leda_polygon ledaP;
    W >> ledaP;
    if (ledaP.size() > 2)
       convert_leda_polygon_to_CGAL(ledaP,p);
}

template <class OutputStream>
void label_vertices(OutputStream& OS, const Polygon_2& poly)
{
   Polygon_2::const_iterator vert_it;

   OS << CGAL::BLACK;
   for (vert_it = poly.vertices_begin(); vert_it != poly.vertices_end(); 
        vert_it++)
   {
      leda_point vertex(CGAL::to_double((*vert_it).x()), 
                        CGAL::to_double((*vert_it).y()));
      leda_string vertex_label("(%.1f, %.1f)", 
                               vertex.xcoord(), vertex.ycoord());
      OS.draw_text(vertex, vertex_label);
   }
}

template <class OutputStream>
void draw_a_polygon(OutputStream& OS, const Polygon_2& poly, CGAL::Color color)
{
   leda_polygon ledaP;
   convert_point_list_to_leda_poly(poly, ledaP);
   OS << color;
   OS.draw_polygon(ledaP);
}

template <class OutputStream>
void draw_poly_list(OutputStream& OS,  const std::list<Polygon_2>& poly_list, 
                    int line_width, CGAL::Color color,  int& count, 
                    bool show_coords)
{
   OS.set_line_width(line_width);
   OS << color;
   std::list<Polygon_2>::const_iterator p_it;

   count = 0;
   for(p_it = poly_list.begin();p_it != poly_list.end(); p_it++)
   {
      Polygon_2::Edge_const_circulator     e_circ = (*p_it).edges_circulator();
      if (e_circ != NULL) 
      {
         Polygon_2::Edge_const_circulator end_circ = e_circ;
         do
         {
           OS << *e_circ;
           OS << (*e_circ).source();
           ++e_circ;
         }
         while (e_circ != end_circ);
      }
      if (show_coords) 
      {
         label_vertices(OS, *p_it);
         OS << color;
      }
      count++;
   }
}

bool read_poly_from_file(Polygon_2& polygon, leda_string file_name)
{
   Point_2  vertex;
   std::ifstream in_file;
   in_file.open(file_name);
   if (in_file) 
   {
       polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
       int num_pieces;
       in_file >> num_pieces;
       in_file >> polygon;
       return true;
   }
   return false;
}

bool save_partition_p_list_to_file(const std::list<Polygon_2>& partition_p_list,
                                   leda_string& current_file_name, 
                                   leda_string file_suffix)
{
   std::ofstream out_file;

   leda_string file_name = current_file_name(current_file_name.pos("/")+1,
                                             current_file_name.length()-1);
   file_name = "output/" + file_name + file_suffix;
   out_file.open(file_name);

   if (!out_file) return false;

   out_file << partition_p_list.size() << std::endl;
   std::list<Polygon_2>::const_iterator it;

   for (it = partition_p_list.begin(); it != partition_p_list.end(); it++) 
   {
      out_file << *it << std::endl;
   }
   out_file.close();
   return true;
}


bool save_poly_to_file(const Polygon_2& polygon, leda_string file_name)
{
   std::ofstream out_file;
   out_file.open(file_name);
   if (out_file)
   {
      out_file << 1 << std::endl;
      out_file << polygon << std::endl;
      out_file.close();
      return true;
   }
   return false;
}

