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
// file          : demo/GeoWin/CGAL_constrained_triang_2.C
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
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/geowin_support.h>
#include <CGAL/leda_rational.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

//#if !defined(_MSC_VER)
//typedef leda_rational coord_type;
//#else
typedef double coord_type;
//#endif

typedef CGAL::Cartesian<coord_type>                             K;
typedef K::Point_2                                              Point;
typedef K::Segment_2                                            Segment;

typedef CGAL::Triangulation_euclidean_traits_2<K>               Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt>                   Vb;
typedef CGAL::Triangulation_face_base_2<Gt>                     Fb;
typedef CGAL::Constrained_triangulation_face_base_2<Gt>         CFb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,CFb> Tds;
typedef CGAL::Constrained_triangulation_2<Gt,Tds>               Constr_triangulation_2;
typedef Constr_triangulation_2::Constraint                      Constraint;
typedef Constr_triangulation_2::Edge                            Edge;
typedef Constr_triangulation_2::Edge_iterator                   Edge_iterator;


class constr_tria : public geowin_update<std::list<Segment>, std::list<Segment> >
{
public:
 void update(const std::list<Segment>& Li, std::list<Segment>& Sl)
 {

  Sl.clear();      
  std::vector<Constraint> SC;
  std::list<Segment>::const_iterator it;
  
  it=Li.begin();

  for(; it != Li.end(); ++it) {
       Segment sakt= *it;
       Point pa = sakt.source(), pb = sakt.target();
       std::pair<Point,Point> pact(pa,pb);
       SC.push_back(pact);
  }

  Constr_triangulation_2 ct(SC.begin(), SC.end() );   
     
  Edge_iterator eit = ct.edges_begin();
  Edge_iterator beyond = ct.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;         
       Sl.push_back(ct.segment(eact));               
       ++eit;  
  } 
 }
};

int main()
{
  geowin_init_default_type((std::list<Segment>*)0, leda_string("CGALSegmentList"));

  std::list<Segment> SL;

  GeoWin GW("CGAL - Constrained Triangulation demo");
  GW.message("The segments of your input should not intersect.");
 
  geo_scene my_scene= GW.new_scene(SL); 
  GW.set_color(my_scene, leda_blue); 
  GW.set_active_line_width(my_scene, 2);

  constr_tria CT;
  geo_scene res2 = GW.new_scene(CT, my_scene, leda_string("Constrained Triangulation"));
  GW.set_color(res2, leda_red);
  GW.set_visible(res2, true);
  
  GW.edit(my_scene);
  
  return 0;  
}

#endif
