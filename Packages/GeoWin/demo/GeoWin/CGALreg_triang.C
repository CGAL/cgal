#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 400)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.0 or higher installed!\n";
 std::cout << "A LEDA version >= 4.0 is required to run GeoWin!\n";
 return 0;
}
#else 

#include <CGAL/Cartesian.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

typedef CGAL::Cartesian<double>  Rep;
typedef CGAL::Circle_2<Rep>       Circle;
typedef CGAL::Point_2<Rep>        Point;

typedef double W;

typedef CGAL::Regular_triangulation_euclidean_traits_2<Rep,W> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Regular_triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt,Tds>  Triangulation_2;
typedef CGAL::Regular_triangulation_2<Gt,Tds> Regular_triangulation_2;

typedef Gt::Weighted_point Weighted_point;

typedef Regular_triangulation_2::Edge Edge;
typedef Regular_triangulation_2::Locate_type Locate_type;
typedef Regular_triangulation_2::Edge_iterator  Edge_iterator;

#include <CGAL/geowin_support.h>

class geo_reg_triang : public geowin_update<std::list<CGALCircle>, std::list<CGALSegment> >
{
public:
 void update(const CGALCirclelist& L, CGALSegmentlist& Sl)
 {
  Regular_triangulation_2 tr;    
  Sl.clear();      
                   
  std::list<CGALCircle>::const_iterator it;
  it= L.begin();
  CGALCircle cakt;
 
  for (; it != L.end() ; ++it) { 
    cakt= *it; 
    //std::cout << Weighted_point(cakt.center(),cakt.squared_radius()) << "\n";
    tr.insert(Weighted_point(cakt.center(),cakt.squared_radius())); 
  }

  Edge_iterator eit = tr.edges_begin();
  Edge_iterator beyond = tr.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;         
       Sl.push_back(tr.segment(eact));               
       ++eit;  
  }         
 }
};

int main()
{
  geowin_init_default_type((CGALCirclelist*)0, leda_string("CGALCircleList"));

  CGALCirclelist L;

  GeoWin GW("CGAL - Regular Triangulation demo");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_active_line_width(my_scene, 2);
  GW.set_color(my_scene, leda_pink);

  geo_reg_triang triangulate;
  geo_scene res1 = GW.new_scene(triangulate ,my_scene , "Triangulation");
  GW.set_color(res1, leda_blue);
  GW.set_line_width(res1, 2);
  GW.set_all_visible(true);
  
  GW.edit(my_scene);
  
  return 0;  
}

#endif
