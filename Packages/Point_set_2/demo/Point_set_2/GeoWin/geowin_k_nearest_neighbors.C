#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 400)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.0 or higher installed!\n";
 std::cout << "A LEDA version >= 4.0 is required !\n";
 return 0;
}
#else 


#include <CGAL/geowin_support.h>
#include <CGAL/Point_set_2.h>

typedef CGAL::Cartesian<double>          REP;
typedef CGAL::point_set_traits_2<REP>      TRAITS;
typedef CGAL::Point_set_2<TRAITS>::Edge    Edge;
typedef CGAL::Point_set_2<TRAITS>::Vertex  Vertex;

CGAL::Point_set_2<TRAITS> PST;
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

class mst : public geowin_update<std::list<CGAL::Point_2<REP> >, std::list<CGAL::Segment_2<REP> > >
{
public:
 void update(const std::list<CGAL::Point_2<REP> >& Lin, std::list<CGAL::Segment_2<REP> >& Lout)
 {
  Lout.clear();
  std::list<Edge> output;  
  std::list<Edge>::const_iterator pit;
  
  PST.minimum_spanning_tree( std::back_inserter(output));  
  
  for (pit=output.begin(); pit != output.end(); pit++){
    Lout.push_back(CGALSegment(PST.seg(*pit)));
  } 
 }
};

class nearest_neighbors : public geowin_update<std::list<CGAL::Point_2<REP> >, std::list<CGAL::Circle_2<REP> > >
{
public:
 void update(const std::list<CGAL::Point_2<REP> >& Lin, std::list<CGAL::Circle_2<REP> >& Lout)
 {
  Lout.clear();
  std::list<CGAL::Point_2<REP> >::const_iterator it = Lin.begin();
  std::list<Vertex> output;  
  std::list<Vertex>::const_iterator pit;

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
  
  mst MS;
  geo_scene sc5 = gw.new_scene(MS, sc1, leda_string("Minimum spanning tree"));
  gw.set_line_width(sc5,2);
 
  gw.set_all_visible(true);
  gw.add_dependence(sc1,sc4);
   
  gw.edit(sc1);

  return 0;
}
#endif
