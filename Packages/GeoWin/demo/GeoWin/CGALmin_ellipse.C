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
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/geowin_support.h>

typedef  CGAL::Cartesian<double>           R;
typedef  CGAL::Point_2<R>                    Point;
typedef  CGAL::Min_ellipse_2_traits_2< R >   Traits;
typedef  CGAL::Min_ellipse_2< Traits >       Min_ellipse;

class geo_ellipse : public geowin_redraw, public geowin_update<CGALPointlist, CGALPointlist >
{
public:
  virtual ~geo_ellipse() {}

  Min_ellipse min_ell;

  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {  W.set_color(c1); W << min_ell.ellipse(); }  

  virtual void update(const CGALPointlist& L, CGALPointlist&)
  {  min_ell.clear(); min_ell.insert(L.begin(),L.end()); }
};

int main()
{
  geowin_init_default_type((CGALPointlist*)0, leda_string("CGALPointList"));

  CGALPointlist L;
  GeoWin GW("CGAL - Optimisation demo - minimal enclosing ellipse");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_point_style(my_scene, leda_circle_point);
  
  geo_ellipse EL;
  geo_scene result  = GW.new_scene(EL, EL , my_scene, leda_string("Minimal enclosing ellipse"));
  GW.set_color(result, leda_red);
  GW.set_visible(result,true);

  GW.edit();
  return 0;
}

#endif
