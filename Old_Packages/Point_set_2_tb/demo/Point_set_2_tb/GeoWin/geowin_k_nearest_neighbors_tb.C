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
int k;

class construct_pointset : public geowin_update<std::list<CGAL::Point_2<REP> >,std::list<CGAL::Segment_2<REP> > >
{
public:
 void update(const std::list<CGAL::Point_2<REP> >& Lin, std::list<CGAL::Segment_2<REP> >& Lout)
 {
  PST.init(Lin.begin(),Lin.end());
  Lout.clear();
  PST.segments(std::back_inserter(Lout));
 }
};


class nearest_neighbors : public geowin_update<std::list<CGAL::Point_2<REP> >, std::list<CGAL::Circle_2<REP> > >
{
public:
 void update(const std::list<CGAL::Point_2<REP> >& Lin, std::list<CGAL::Circle_2<REP> >& Lout)
 {
  Lout.clear();
  std::list<CGAL::Point_2<REP> >::const_iterator it = Lin.begin();
  std::list<Vertex_handle> output;  
  std::list<Vertex_handle>::const_iterator pit;

  for (; it != Lin.end(); it++){
    PST.nearest_neighbors(*it, k, std::back_inserter(output));
  } 
  for (pit=output.begin(); pit != output.end(); pit++){
    Lout.push_back(CGALCircle(PST.pos(*pit),2.0));
  } 
 }
};

int main()
{
  geowin_init_default_type((std::list<CGAL::Point_2<REP> >*)0, leda_string("CGALPointlist"));
    
  std::cout << "Find the k nearest neighbors of every point in scene 2.\n";
  std::cout << "k:"; std::cin >> k;
  
  GeoWin gw;

  std::list<CGAL::Point_2<REP> > Lp;
  geo_scene sc1 = gw.new_scene(Lp);
  
  std::list<CGAL::Point_2<REP> > Lother;
  geo_scene sc2 = gw.new_scene(Lother);
  gw.set_color(sc2,leda_blue);
  
  construct_pointset CP;
  geo_scene sc3 = gw.new_scene(CP, sc1, leda_string("2d point set"));
  gw.set_color(sc3,leda_blue);
  
  nearest_neighbors NN;
  geo_scene sc4 = gw.new_scene(NN, sc2, leda_string("k nearest neighbors"));
  gw.set_fill_color(sc4,leda_red);
  gw.set_color(sc4,leda_red);
  gw.set_point_style(sc4,leda_circle_point);
  gw.set_line_width(sc4,4);
 
  gw.set_all_visible(true);
  gw.add_dependence(sc1,sc4);
   
  gw.edit(sc1);

  return 0;
}

