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
#include <CGAL/leda_integer.h>

#if defined GEOWIN_SUPPORT_NO_TEMPLATES
#include "geowin_workaround1.h"
#endif

typedef CGAL::Cartesian<leda_integer>   My_REP;
typedef CGAL::Point_2< My_REP >     My_CGALPoint;
typedef std::list<My_CGALPoint>     My_CGALPointlist;
typedef CGAL::Segment_2< My_REP >   My_CGALSegment;
typedef std::list<My_CGALSegment>   My_CGALSegmentlist;


class geo_hull : public geowin_update<std::list<My_CGALPoint>, std::list<My_CGALSegment> >
{
public:
 void update(const My_CGALPointlist& L, My_CGALSegmentlist& Sl)
 {
  Sl.clear();
  My_CGALPointlist out;

  CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(out));   

  // building the segment list ...
  if( out.size() > 1 ) {
    My_CGALPoint pakt,prev,pstart;

    std::list<My_CGALPoint>::const_iterator it;
    it=out.begin();
    prev= *it; pstart=prev;
    it++;

    for(; it != out.end(); ++it) {
       pakt= *it;
       Sl.push_back(My_CGALSegment(prev,pakt));
       prev=pakt;
    }

    Sl.push_back(My_CGALSegment(pakt,pstart));
  }
 }
};

int main()
{
  geowin_init_default_type((My_CGALPointlist*)0, leda_string("My_CGALPointList"));
 
  My_CGALPointlist L;
  My_CGALPointlist out;

  GeoWin GW("CGALTEST - convex hull");

  geo_hull update_obj;
  geo_scene my_scene= GW.new_scene(L);  
  geo_scene result  = GW.new_scene(update_obj, my_scene, leda_string("Convex Hull"));
   
  GW.set_visible(result,true);
  GW.edit(my_scene);
  
  return 0;  
}

#endif
