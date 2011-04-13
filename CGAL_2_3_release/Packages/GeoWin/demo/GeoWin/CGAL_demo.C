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

int main()
{

  geowin_init_default_type((CGALPointlist*)0, leda_string("CGALPointList"));
  geowin_init_default_type((CGALSegmentlist*)0, leda_string("CGALSegmentList"));
  geowin_init_default_type((CGALCirclelist*)0, leda_string("CGALCircleList"));
  geowin_init_default_type((CGALLinelist*)0, leda_string("CGALLineList"));
  geowin_init_default_type((CGALRaylist*)0, leda_string("CGALRayList"));
  geowin_init_default_type((CGALTrianglelist*)0, leda_string("CGALTriangleList"));
  geowin_init_default_type((CGALRectanglelist*)0, leda_string("CGALRectangleList"));

  geowin_init_default_type((CGALPolygonlist*)0, leda_string("CGALPolygonList"));

  geowin_init_default_type((CGALPoint_3_list*)0, leda_string("CGALPoint_3_List")); 

  GeoWin GW("GeoWin Demo using CGAL objects");
  GW.edit();
  return 0;  
}

#endif
