#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 420)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.2 or higher installed!\n";
 std::cout << "A LEDA version >= 4.2 is required !\n";
 return 0;
}
#else 

#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits      K;
typedef K::Point_2                        Point;
typedef K::Segment_2                      Segment;

class geo_hull : public geowin_update<std::list<Point>, std::list<Segment> >
{
public:
 void update(const std::list<Point>& L, std::list<Segment>& Sl)
 {
  Sl.clear();
  std::list<Point> out;
  
  K traits;
  CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(out), traits);   

  if( out.size() > 1 ) {
    Point pakt,prev,pstart;

    std::list<Point>::const_iterator it=out.begin();
    prev= *it; pstart=prev;
    it++;

    for(; it != out.end(); ++it) {
       pakt= *it;
       Sl.push_back(Segment(prev,pakt));
       prev=pakt;
    }
    Sl.push_back(Segment(pakt,pstart));
  }
 }
 
};

int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("LEDA-rat_point"));
 
  std::list<Point> L;

  GeoWin GW("2d convex hull");
  geo_hull update_obj;
  geo_scene my_scene= GW.new_scene(L);  
  geo_scene result  = GW.new_scene(update_obj,my_scene,leda_string("Convex Hull")); 
  GW.set_visible(result,true);
 
  GW.edit(my_scene);
  
  return 0;  
}

#endif
