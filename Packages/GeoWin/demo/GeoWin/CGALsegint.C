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
#include <CGAL/Segment_2_Segment_2_intersection.h>
#include <CGAL/geowin_support.h>

class geo_segint : public geowin_update<std::list<CGALSegment>,std::list<CGALPoint> >
{
public:
 void update(const CGALSegmentlist& L, CGALPointlist& Sl)
 {
  Sl.clear();
  if (L.size() <2) return;
  int cnt=0;
  CGALSegment s1,s2;
  CGALSegmentlist::const_iterator it1= L.begin(), it2;
  CGAL::Object result;
  CGALPoint pt;

  for(;it1 != L.end(); ++it1){
   for(it2=L.begin() ;it2 != L.end(); ++it2){
    s1= *it1; s2= *it2;
    result=CGAL::intersection(s1,s2);
    if (CGAL::assign(pt,result) ) { Sl.push_back(pt); cnt++; }
   }
  }
 }
};

int main()
{
  geowin_init_default_type((CGALSegmentlist*)0, leda_string("CGALSegmentList"));
 
  CGALSegmentlist L;

  GeoWin GW("CGAL - Segment Intersection");
 
  geo_scene my_scene= GW.new_scene(L);  

  geo_segint inter;
  geo_scene res  = GW.new_scene(inter,my_scene,leda_string("Segment intersection"));
  GW.set_color( res, leda_red );
  GW.set_point_style( res, leda_disc_point );
  GW.set_fill_color( res, leda_red);
  GW.set_visible(res,true);
 
  GW.edit(my_scene);
  
  return 0;  
}

#endif
