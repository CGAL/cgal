#include <CGAL/geowin_support.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>

#include <CGAL/Point_set_2_tb.h>

typedef CGAL::Cartesian<double>          REP;
typedef CGAL::Triangulation_euclidean_traits_2<REP> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;

typedef CGAL::point_set_traits_2<REP>      TRAITS;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge    Edge;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex  Vertex;


CGAL::Point_set_2_tb<TRAITS,Gt,Tds> PST;
  
class construct_pointset : public geowin_update<std::list<CGALPoint>, std::list<CGALSegment> >
{
public:
 void update(const std::list<CGALPoint>& Lin, std::list<CGALSegment>& Lout)
 {
  PST.init(Lin.begin(),Lin.end());
  Lout.clear();
  PST.segments(std::back_inserter(Lout));
 }
};


class circle_range_search : public geowin_update<std::list<CGALCircle>, std::list<CGALCircle> >
{
 void update(const std::list<CGALCircle>& Lin, std::list<CGALCircle>& Lout)
 {
  Lout.clear();
  std::list<CGALCircle>::const_iterator it = Lin.begin();
  std::list<Vertex_handle> output;  
  std::list<Vertex_handle>::const_iterator pit;

  for (; it != Lin.end(); it++){
    PST.range_search(*it, std::back_inserter(output));
  } 
  for (pit=output.begin(); pit != output.end(); pit++){
    Lout.push_back(CGALCircle(PST.pos(*pit),2.0));
  } 
 }
};

class triangle_range_search : public geowin_update<std::list<CGALTriangle>, std::list<CGALCircle> >
{
public:
 void update(const std::list<CGALTriangle>& Lin, std::list<CGALCircle>& Lout)
 {
  Lout.clear();
  std::list<CGALTriangle>::const_iterator it = Lin.begin();
  std::list<Vertex_handle> output;  
  std::list<Vertex_handle>::const_iterator pit;

  for (; it != Lin.end(); it++){
    PST.range_search((*it).vertex(0),(*it).vertex(1),(*it).vertex(2), std::back_inserter(output));
  } 
  for (pit=output.begin(); pit != output.end(); pit++){
    Lout.push_back(CGALCircle(PST.pos(*pit),1.0));
  } 
 }
};


int main()
{
  geowin_init_default_type((CGALPointlist*)0, leda_string("CGALPointlist"));
  geowin_init_default_type((CGALCirclelist*)0, leda_string("CGALCirclelist"));
  geowin_init_default_type((CGALTrianglelist*)0, leda_string("CGALTrianglelist"));
    
  GeoWin gw;

  std::list<CGALPoint> Lp;
  geo_scene sc1 = gw.new_scene(Lp);
 
  std::list<CGALCircle> Lc;
  geo_scene sc_circ = gw.new_scene(Lc);
  gw.set_color(sc_circ, leda_blue);
  gw.set_fill_color(sc_circ, leda_invisible);
  
  std::list<CGALTriangle> Lt;
  geo_scene sc_triang = gw.new_scene(Lt);
  gw.set_color(sc_triang, leda_blue2);
  gw.set_line_width(sc_triang,3);
  gw.set_active_line_width(sc_triang,3); 
  gw.set_fill_color(sc_triang, leda_invisible);
  
  construct_pointset CP;
  gw.new_scene(CP, sc1, leda_string("2d point set"));
  
  circle_range_search CRS;
  geo_scene sc3 = gw.new_scene(CRS, sc_circ, leda_string("circular range search"));
  gw.set_color(sc3,leda_red);
  gw.set_point_style(sc3, leda_circle_point);
  gw.set_line_width(sc3,4);
  
  triangle_range_search TRS;
  geo_scene sc4 = gw.new_scene(TRS, sc_triang, leda_string("triangular range search"));
  gw.set_color(sc4,leda_green);
  gw.set_point_style(sc4, leda_circle_point);  
  gw.set_line_width(sc4,4);
  
  gw.set_all_visible(true);
  gw.add_dependence(sc1,sc3);
  gw.add_dependence(sc1,sc4);

  gw.edit(sc1);

  return 1;
}

