// ======================================================================
//
// Copyright (c) 1999 The GALIA Consortium
//
// This software and related documentation is part of the
// Computational Geometry Algorithms Library (CGAL).
//
// Every use of CGAL requires a license. Licenses come in three kinds:
//
// - For academic research and teaching purposes, permission to use and
//   copy the software and its documentation is hereby granted free of  
//   charge, provided that
//   (1) it is not a component of a commercial product, and
//   (2) this notice appears in all copies of the software and
//       related documentation.
// - Development licenses grant access to the source code of the library 
//   to develop programs. These programs may be sold to other parties as 
//   executable code. To obtain a development license, please contact
//   the GALIA Consortium (at cgal@cs.uu.nl).
// - Commercialization licenses grant access to the source code and the
//   right to sell development licenses. To obtain a commercialization 
//   license, please contact the GALIA Consortium (at cgal@cs.uu.nl).
//
// This software and documentation is provided "as-is" and without
// warranty of any kind. In no event shall the CGAL Consortium be
// liable for any damage of any kind.
//
// The GALIA Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Trier University
// (Germany), Max-Planck-Institute Saarbrucken (Germany),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : demo/GeoWin/CGAL_polygon_2.C
//
// ======================================================================

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

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::Cartesian<double>                      K;
typedef K::Point_2                                   Point;
typedef CGAL::Polygon_traits_2<K>                    PTraits;
typedef CGAL::Polygon_2<PTraits,std::list<Point> >   Polygon;


class geo_poly_loc : public geowin_redraw, public geowin_update<std::list<Polygon>, std::list<Point> >
{
public:
  geo_scene input_points; 
  std::list<Point> in; 
  std::list<Point> out;
  std::list<Point> bound;

  virtual ~geo_poly_loc() {}

  virtual bool write_postscript(ps_file& PS,leda_color c1,leda_color c2) 
  { 
   std::list<Point>::const_iterator rt;

   PS.set_color(c1);
   for (rt=in.begin(); rt != in.end(); rt++) PS.draw_disc((*rt).x(),(*rt).y(),5);
   PS.set_color(c2);
   for (rt=out.begin(); rt != out.end(); rt++) PS.draw_disc((*rt).x(),(*rt).y(),5);
   for (rt=bound.begin(); rt != bound.end(); rt++) PS.draw_disc((*rt).x(),(*rt).y(),5);
   
    return true;
  }

  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {  
   std::list<Point>::const_iterator rt;

   W.set_color(c1);
   for (rt=in.begin(); rt != in.end(); rt++) W.draw_disc((*rt).x(),(*rt).y(),5);
   W.set_color(c2);
   for (rt=out.begin(); rt != out.end(); rt++) W.draw_disc((*rt).x(),(*rt).y(),5);
   for (rt=bound.begin(); rt != bound.end(); rt++) W.draw_disc((*rt).x(),(*rt).y(),5);
  }

  virtual void update(const std::list<Polygon>& L, std::list<Point>&)
  { 
    if (L.empty()) return;
    Polygon polygon = (*(L.begin()));
    in.clear(); out.clear(); bound.clear();

    GeoWin* gw= get_geowin(input_points);
    std::list<Point> LST;
    gw->get_objects(input_points,LST);
 
    std::list<Point>::const_iterator it;
    Point lakt;

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
  geowin_init_default_type((std::list<Point>*)0, leda_string("CGALPointList"));
  geowin_init_default_type((std::list<Polygon>*)0, leda_string("CGALPolygonList"));
 
  std::list<Polygon>  L;
  std::list<Point>    L2;

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
