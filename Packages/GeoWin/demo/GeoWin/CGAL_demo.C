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
// file          : demo/GeoWin/CGAL_demo.C
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
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::Cartesian<double>                      K;
typedef K::Point_2                                   Point;
typedef K::Segment_2                                 Segment;
typedef K::Circle_2                                  Circle;
typedef K::Line_2                                    Line;
typedef K::Ray_2                                     Ray;
typedef K::Triangle_2                                Triangle;
typedef K::Iso_rectangle_2                           Iso_rectangle;
typedef K::Point_3                                   Point_3;
typedef CGAL::Polygon_traits_2<K>                    PTraits;
typedef CGAL::Polygon_2<PTraits,std::list<Point> >   Polygon;

int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("CGALPointList"));
  geowin_init_default_type((std::list<Segment>*)0, leda_string("CGALSegmentList"));
  geowin_init_default_type((std::list<Circle>*)0, leda_string("CGALCircleList"));
  geowin_init_default_type((std::list<Line>*)0, leda_string("CGALLineList"));
  geowin_init_default_type((std::list<Ray>*)0, leda_string("CGALRayList"));
  geowin_init_default_type((std::list<Triangle>*)0, leda_string("CGALTriangleList"));
  geowin_init_default_type((std::list<Iso_rectangle>*)0, leda_string("CGALRectangleList"));
  geowin_init_default_type((std::list<Polygon>*)0, leda_string("CGALPolygonList"));
  geowin_init_default_type((std::list<Point_3>*)0, leda_string("CGALPoint_3_List")); 

  GeoWin GW("GeoWin Demo using CGAL objects");
  GW.edit();
  return 0;  
}

#endif
