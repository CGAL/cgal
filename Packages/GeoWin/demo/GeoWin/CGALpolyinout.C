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
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/geowin_support.h>

class geo_poly_loc : public geowin_redraw, public geowin_update<CGALPolygonlist, CGALPointlist >
{
public:
  geo_scene input_points; 
  std::list<CGALPoint> in; 
  std::list<CGALPoint> out;
  std::list<CGALPoint> bound;

  virtual ~geo_poly_loc() {}

  virtual bool write_postscript(ps_file& PS,leda_color c1,leda_color c2) 
  { 
   std::list<CGALPoint>::const_iterator rt;

   PS.set_color(c1);
   for (rt=in.begin(); rt != in.end(); rt++) PS.draw_disc((*rt).x(),(*rt).y(),5);
   PS.set_color(c2);
   for (rt=out.begin(); rt != out.end(); rt++) PS.draw_disc((*rt).x(),(*rt).y(),5);
   for (rt=bound.begin(); rt != bound.end(); rt++) PS.draw_disc((*rt).x(),(*rt).y(),5);
   
    return true;
  }

  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {  
   std::list<CGALPoint>::const_iterator rt;

   W.set_color(c1);
   for (rt=in.begin(); rt != in.end(); rt++) W.draw_disc((*rt).x(),(*rt).y(),5);
   W.set_color(c2);
   for (rt=out.begin(); rt != out.end(); rt++) W.draw_disc((*rt).x(),(*rt).y(),5);
   for (rt=bound.begin(); rt != bound.end(); rt++) W.draw_disc((*rt).x(),(*rt).y(),5);
  }

  virtual void update(const CGALPolygonlist& L, CGALPointlist&)
  { 
    if (L.empty()) return;
    CGALPolygon polygon = (*(L.begin()));
    in.clear(); out.clear(); bound.clear();

    GeoWin* gw= get_geowin(input_points);
    CGALPointlist LST;
    gw->get_objects(input_points,LST);
 
    std::list<CGALPoint>::const_iterator it;
    CGALPoint lakt;

    for(it=LST.begin(); it != LST.end(); ++it) { 
      lakt=*it;
      CGAL::Bounded_side bside   = polygon.bounded_side(lakt);
      switch (bside) {
        case CGAL::ON_BOUNDED_SIDE:
          in.push_back(lakt); break;

        case CGAL::ON_BOUNDARY:
          bound.push_back(lakt); break;

        case CGAL::ON_UNBOUNDED_SIDE:
          out.push_back(lakt); break;
      }
    }

  }
};


int main()
{
  geowin_init_default_type((CGALPointlist*)0, leda_string("CGALPointList"));
  geowin_init_default_type((CGALPolygonlist*)0, leda_string("CGALPolygonList"));
 
  CGALPolygonlist L;
  CGALPointlist L2;

  GeoWin GW("CGAL - inside polygon test");
 
  geo_scene poly_scene= GW.new_scene(L);  
  GW.set_color(poly_scene, leda_red );
  GW.set_line_width(poly_scene, 3);
  GW.set_active_line_width(poly_scene, 3);
  GW.set_fill_color(poly_scene, leda_green2);

  geo_scene point_scene= GW.new_scene(L2);  
  GW.set_color( point_scene, leda_red );
 
  geo_poly_loc RU;
  RU.input_points = point_scene;

  geo_scene RES = GW.new_scene( RU, RU, poly_scene, leda_string("inside/outside")); 
  GW.set_color( RES, leda_green );

  GW.set_all_visible(true);
  GW.add_dependence(point_scene,RES);
 
  GW.edit(poly_scene);
  
  return 0;  
}

#endif
