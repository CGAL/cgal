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
// file          : demo/GeoWin/CGAL_convex_hull_2.C
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
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/geowin_support.h>

class geo_hull : public geowin_update<std::list<CGALPoint>, std::list<CGALSegment> >
{
public:
 void update(const CGALPointlist& L, CGALSegmentlist& Sl)
 {
  Sl.clear();
  CGALPointlist out;
  CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(out));   

  if( out.size() > 1 ) {
    CGALPoint pakt,prev,pstart;

    std::list<CGALPoint>::const_iterator it;
    it=out.begin();
    prev= *it; pstart=prev;
    it++;

    for(; it != out.end(); ++it) {
       pakt= *it;
       Sl.push_back(CGALSegment(prev,pakt));
       prev=pakt;
    }
    Sl.push_back(CGALSegment(pakt,pstart));
  }
 }
 
};

int main()
{
  geowin_init_default_type((CGALPointlist*)0, leda_string("CGALPointList"));
 
  CGALPointlist L;
  CGALPointlist out;

  GeoWin GW("CGALTEST - convex hull");

  geo_hull update_obj;
  geo_scene my_scene= GW.new_scene(L);  
  geo_scene result  = GW.new_scene(update_obj,my_scene,leda_string("Convex Hull")); 
  GW.set_visible(result,true);
 
  GW.edit(my_scene);
  
  return 0;  
}

#endif
