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
// file          : demo/GeoWin/CGAL_ler_2.C
//
// ======================================================================
// Largest empty rectangle demo

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
#include <CGAL/Largest_empty_iso_rectangle_2.h>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef double                                                  NT;
typedef CGAL::Cartesian<NT>                                     K;
typedef CGAL::Largest_empty_iso_rectangle_2<K>                  LER;
typedef K::Iso_rectangle_2                                      Iso_rectangle;
typedef K::Point_2                                              Point;

std::list<Iso_rectangle> bbox; // stores the bounding box ...

class geo_ler : public geowin_update<std::list<Point>, std::list<Iso_rectangle> >
{
public:
 void update(const std::list<Point>& L, std::list<Iso_rectangle>& R)
 {
  R.clear(); 
  if (bbox.size() == 0) return;
  LER largest_empty_rec(* bbox.begin());
  largest_empty_rec.insert(L.begin(),L.end());   
  Iso_rectangle result = largest_empty_rec.get_largest_empty_iso_rectangle();
  R.push_back(result);
 }
};

int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("CGALPointList"));
  geowin_init_default_type((std::list<Iso_rectangle>*)0, leda_string("CGALIsorectangleList"));
  
  std::list<Point> L;
  GeoWin GW("2d Largest empty rectangle in a set of points");
  GW.message("You have to input a bounding box in the isorectangle scene");

  geo_ler update_obj;
  geo_scene my_scene= GW.new_scene(L);  
  geo_scene rect_scene= GW.new_scene(bbox);
  GW.set_fill_color(rect_scene, leda_invisible);
  
  geo_scene result  = GW.new_scene(update_obj,my_scene,leda_string("LER")); 
  GW.set_color(result, leda_red);
  GW.set_fill_color(result, leda_red);
  GW.set_line_width(result, 2);
   
  GW.add_dependence(rect_scene, result);
  GW.set_all_visible(true);
  GW.edit(my_scene); 
  return 0;  
}
#endif
