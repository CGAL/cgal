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
// file          : demo/GeoWin/CGAL_min_ellipse_2.C
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
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef  CGAL::Cartesian<double>             K;
typedef  K::Point_2                          Point;
typedef  CGAL::Min_ellipse_2_traits_2<K>     Traits;
typedef  CGAL::Min_ellipse_2< Traits >       Min_ellipse;

class geo_ellipse : public geowin_redraw, public geowin_update<std::list<Point>,std::list<Point> >
{
public:
  virtual ~geo_ellipse() {}

  Min_ellipse min_ell;

  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {  W.set_color(c1); W << min_ell.ellipse(); }  

  virtual void update(const std::list<Point>& L, std::list<Point>&)
  {  min_ell.clear(); min_ell.insert(L.begin(),L.end()); }
};

int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("CGALPointList"));

  std::list<Point> L;
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
