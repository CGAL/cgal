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
// file          : demo/GeoWin/CGAL_segment_intersection_2.C
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
#include <CGAL/Segment_2_Segment_2_intersection.h>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::Cartesian<double>   K;
typedef K::Segment_2              Segment;
typedef K::Point_2                Point;

class geo_segint : public geowin_update<std::list<Segment>,std::list<Point> >
{
public:
 void update(const std::list<Segment>& L, std::list<Point>& Sl)
 {
  Sl.clear();
  if (L.size() <2) return;
  int cnt=0;
  Segment s1,s2;
  std::list<Segment>::const_iterator it1= L.begin(), it2;
  CGAL::Object result;
  Point pt;

  for(;it1 != L.end(); ++it1){
   for(it2=L.begin() ;it2 != L.end(); ++it2){
    s1= *it1; s2= *it2;
    result=CGAL::intersection(s1,s2);
    if (CGAL::assign(pt,result) ) { Sl.push_back(pt); cnt++; }
   }
  }
 }
};

int main()
{
  geowin_init_default_type((std::list<Segment>*)0, leda_string("CGALSegmentList"));
 
  std::list<Segment> L;

  GeoWin GW("CGAL - Segment Intersection");
 
  geo_scene my_scene= GW.new_scene(L);  

  geo_segint inter;
  geo_scene res  = GW.new_scene(inter,my_scene,leda_string("Segment intersection"));
  GW.set_color( res, leda_red );
  GW.set_point_style( res, leda_disc_point );
  GW.set_fill_color( res, leda_red);
  GW.set_visible(res,true);
 
  GW.edit(my_scene);
  
  return 0;  
}

#endif
