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
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/geowin_support.h>

class geo_hull : public geowin_update<std::list<CGALPoint>, std::list<CGALSegment> >
{
public:
 void update(const CGALPointlist& L, CGALSegmentlist& Sl)
 {
  Sl.clear();
  CGALPointlist out;
  CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(out));   

  if( out.size() > 1 ) {
    CGALPoint pakt,prev,pstart;

    std::list<CGALPoint>::const_iterator it;
    it=out.begin();
    prev= *it; pstart=prev;
    it++;

    for(; it != out.end(); ++it) {
       pakt= *it;
       Sl.push_back(CGALSegment(prev,pakt));
       prev=pakt;
    }
    Sl.push_back(CGALSegment(pakt,pstart));
  }
 }
 
};

int main()
{
  geowin_init_default_type((CGALPointlist*)0, leda_string("CGALPointList"));
 
  CGALPointlist L;
  CGALPointlist out;

  GeoWin GW("CGALTEST - convex hull");

  geo_hull update_obj;
  geo_scene my_scene= GW.new_scene(L);  
  geo_scene result  = GW.new_scene(update_obj,my_scene,leda_string("Convex Hull")); 
  GW.set_visible(result,true);
 
  GW.edit(my_scene);
  
  return 0;  
}

#endif
