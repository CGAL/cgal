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

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <list>
#include <vector>
#include <CGAL/geowin_support.h>

typedef CGAL::Cartesian<double> R;
typedef CGAL::Polygon_traits_2<R> Traits;
typedef Traits::Point_2 Point;
typedef std::vector<Point> Container;
typedef CGAL::Polygon_2<Traits,Container> CGALVecPolygon;

class geo_all_furthest_nb : public geowin_update<std::list<CGALPoint>,std::list<CGALSegment> >
{
public:

 void update(const CGALPointlist& L, CGALSegmentlist& Sl)
 {
  Sl.clear();
  CGALVecPolygon P;

  // building convex polygon...
  CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(P));   
  
  if (P.size()<2) return;

  std::list<int> il;

  CGAL::all_furthest_neighbors(P.vertices_begin(), P.vertices_end(), std::back_inserter(il));

  std::list<int>::const_iterator lit= il.begin();
  int z=0;
  
  std::vector<CGALPoint> CT = P.container();

  for(; lit != il.end(); ++lit) {
     CGALPoint p1= CT[z];
     CGALPoint p2= CT[*lit]; 

     Sl.push_back(CGALSegment(p1,p2));
     z++;
  }
 }
};

class conv_hull_seg : public geowin_update<std::list<CGALPoint>,std::list<CGALSegment> >
{
public:

 void update(const CGALPointlist& L, CGALSegmentlist& Sl)
 {
  Sl.clear();
  CGALPointlist out;

   CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(out));   

  // building the segment list ...
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

  GeoWin GW("CGAL - All furthest neighbors");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_point_style(my_scene, leda_disc_point); 
  
  conv_hull_seg CHS;
  geo_scene hull_scene  = GW.new_scene(CHS,my_scene,leda_string("Convex Hull"));
  GW.set_visible(hull_scene,true);

  geo_all_furthest_nb AFN;
  geo_scene result  = GW.new_scene(AFN,my_scene,leda_string("All furthest neighbors"));
  GW.set_color(result,leda_red);
  GW.set_visible(result,true);
 
  GW.message("The furtest neighbors of the vertices of the convex hull are displayed");
  GW.edit(my_scene);
  return 0;  
}

#endif
