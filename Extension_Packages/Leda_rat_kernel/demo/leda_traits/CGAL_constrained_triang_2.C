
#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 430)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.3 or higher installed!\n";
 std::cout << "A LEDA version >= 4.3 is required!\n";
 return 0;
}
#else 

#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>  
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif


typedef CGAL::leda_rat_kernel_traits      K;
typedef K::Point_2                        Point;
typedef K::Segment_2                      Segment;

typedef CGAL::Triangulation_vertex_base_2<K>                    Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>          Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>             TDS;
typedef CGAL::Exact_intersections_tag                           Itag;
typedef CGAL::Constrained_triangulation_2<K,TDS,Itag>           CT;
typedef CGAL::Constrained_triangulation_plus_2<CT>              Constr_triangulation_2;

typedef Constr_triangulation_2::Constraint                      Constraint;
typedef Constr_triangulation_2::Edge                            Edge;
typedef Constr_triangulation_2::Edge_iterator                   Edge_iterator;


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
  geowin_init_default_type((std::list<Segment>*)0, leda_string("LEDA-rat_segment"));

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
