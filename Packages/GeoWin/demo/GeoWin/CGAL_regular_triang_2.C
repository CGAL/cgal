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
// file          : demo/GeoWin/CGAL_regular_triang_2.C
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
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

typedef CGAL::Cartesian<double>   K;
typedef K::Circle_2               Circle;
typedef K::Segment_2              Segment;
typedef K::Point_2                Point;

typedef double W;

typedef CGAL::Regular_triangulation_euclidean_traits_2<K,W>    Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt>                  Vb;
typedef CGAL::Regular_triangulation_face_base_2<Gt>            Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt,Tds>                          Triangulation_2;
typedef CGAL::Regular_triangulation_2<Gt,Tds>                  Regular_triangulation_2;

typedef Gt::Weighted_point                                     Weighted_point;
typedef Regular_triangulation_2::Edge                          Edge;
typedef Regular_triangulation_2::Edge_iterator                 Edge_iterator;

#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

class geo_reg_triang : public geowin_update<std::list<Circle>, std::list<Segment> >
{
public:
 void update(const std::list<Circle>& L, std::list<Segment>& Sl)
 {
  Regular_triangulation_2 tr;    
  Sl.clear();      
                   
  std::list<Circle>::const_iterator it = L.begin();
  Circle cakt;
 
  for (; it != L.end() ; ++it) { 
    cakt= *it; 
    tr.insert(Weighted_point(cakt.center(),cakt.squared_radius())); 
  }

  Edge_iterator eit = tr.edges_begin();
  Edge_iterator beyond = tr.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;         
       Sl.push_back(tr.segment(eact));               
       ++eit;  
  }         
 }
};

int main()
{
  geowin_init_default_type((std::list<Circle>*)0, leda_string("CGALCircleList"));

  std::list<Circle> L;

  GeoWin GW("CGAL - Regular Triangulation demo");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_active_line_width(my_scene, 2);
  GW.set_color(my_scene, leda_pink);

  geo_reg_triang triangulate;
  geo_scene res1 = GW.new_scene(triangulate ,my_scene , "Triangulation");
  GW.set_color(res1, leda_blue);
  GW.set_line_width(res1, 2);
  GW.set_all_visible(true);
  
  GW.edit(my_scene);
  
  return 0;  
}

#endif
