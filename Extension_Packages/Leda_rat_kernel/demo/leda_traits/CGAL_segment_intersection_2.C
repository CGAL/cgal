
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
#include <CGAL/Segment_2_Segment_2_intersection.h>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits      K;
typedef K::Segment_2                      Segment;
typedef K::Point_2                        Point;
typedef K::Intersect_2                    Intersect_2;

class geo_segint : public geowin_update<std::list<Segment>,std::list<Point> >
{
public:
 void update(const std::list<Segment>& L, std::list<Point>& Sl)
 {
  Sl.clear();
  if (L.size() <2) return;
  int cnt=0;
  Segment s1,s2;
  std::list<Segment>::const_iterator it1= L.begin(), it2;
  CGAL::Object result;
  Point pt;
  Intersect_2  inter;

  for(;it1 != L.end(); ++it1){
   for(it2=L.begin() ;it2 != L.end(); ++it2){
    s1= *it1; s2= *it2;
    result= inter(s1,s2);
    if (CGAL::assign(pt,result) ) { Sl.push_back(pt); cnt++; }
   }
  }
 }
};

int main()
{
  geowin_init_default_type((std::list<Segment>*)0, leda_string("LEDA-rat_segment"));
 
  std::list<Segment> L;

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
