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
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/geowin_support.h>

typedef  CGAL::Cartesian<double>       R;
typedef  CGAL::Point_2<R>                Point;
typedef  CGAL::Min_circle_2_traits_2<R>  Traits;
typedef  CGAL::Min_circle_2<Traits>      Min_circle;
typedef  Min_circle::Circle              OptCircle;

class geo_circ : public geowin_update<std::list<CGALPoint>,std::list<CGALCircle> >
{
public:
 void update(const CGALPointlist& L, CGALCirclelist& Cl)
 {
   Cl.clear();
   if (L.size() < 2) return;

   Min_circle  mc1( L.begin(), L.end(), true);
   OptCircle ci= mc1.circle();

   Point ctp=ci.center();
   CGALCircle conv(ctp,ci.squared_radius());
   Cl.push_back(conv); 
 }
};

int main()
{
  geowin_init_default_type((CGALPointlist*)0, leda_string("CGALPointList"));

  CGALPointlist L;
  GeoWin GW("CGAL - Optimisation demo");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_point_style(my_scene, leda_disc_point);
  
  geo_circ min_circ;
  geo_scene result  = GW.new_scene(min_circ ,my_scene , leda_string("Minimal circle"));
  GW.set_color(result,leda_blue);
  GW.set_fill_color(result,leda_invisible);
  GW.set_line_width(result, 3);
  GW.set_visible(result,true);

  GW.edit(my_scene);
  return 0;
}

#endif
