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
// file          : demo/GeoWin/CGAL_min_circle_2.C
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
  GW.message("We compute the minimum enclosing circle of the input point set");

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
