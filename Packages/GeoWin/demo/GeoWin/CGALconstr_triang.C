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
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Constrained_triangulation_2.h>

#include <CGAL/leda_rational.h>

//#if !defined(_MSC_VER)
//typedef leda_rational coord_type;
//#else
typedef double coord_type;
//#endif

typedef CGAL::Cartesian<coord_type>  Rep;

typedef CGAL::Point_2<Rep>  Point;
typedef CGAL::Segment_2<Rep>  Segment;

typedef CGAL::Triangulation_euclidean_traits_2<Rep> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Constrained_triangulation_face_base_2<Gt>        CFb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,CFb> Tds;
typedef CGAL::Constrained_triangulation_2<Gt,Tds> Constr_triangulation_2;
typedef Constr_triangulation_2::Constraint  Constraint;
  
typedef Constr_triangulation_2::Edge Edge;
typedef Constr_triangulation_2::Face_handle  Face_handle;
typedef Constr_triangulation_2::Edge_iterator  Edge_iterator;

#include <CGAL/geowin_support.h>

class constr_tria : public geowin_update<std::list<Segment>, std::list<Segment> >
{
public:
 void update(const std::list<Segment>& Li, std::list<Segment>& Sl)
 {

  Sl.clear();      
  std::vector<Constraint> SC;
  std::list<Segment>::const_iterator it;
  
  it=Li.begin();

  for(; it != Li.end(); ++it) {
       Segment sakt= *it;
       Point pa = sakt.source(), pb = sakt.target();
       std::pair<Point,Point> pact(pa,pb);
       SC.push_back(pact);
  }

  Constr_triangulation_2 ct(SC.begin(), SC.end() );   
     
  Edge_iterator eit = ct.edges_begin();
  Edge_iterator beyond = ct.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;         
       Sl.push_back(ct.segment(eact));               
       ++eit;  
  } 
 }
};

int main()
{
  geowin_init_default_type((std::list<Segment>*)0, leda_string("CGALSegmentList"));

  std::list<Segment> SL;

  GeoWin GW("CGAL - Constrained Triangulation demo");
 
  geo_scene my_scene= GW.new_scene(SL); 
  GW.set_color(my_scene, leda_blue); 
  GW.set_active_line_width(my_scene, 2);

  constr_tria CT;
  geo_scene res2 = GW.new_scene(CT, my_scene, leda_string("Constrained Triangulation"));
  GW.set_color(res2, leda_red);
  GW.set_visible(res2, true);
  
  GW.edit(my_scene);
  
  return 0;  
}

#endif
