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
// file          : demo/GeoWin/CGAL_partition_2.C
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
#include <CGAL/Polygon_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/geowin_support.h>

#include <CGAL/leda_rational.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

//typedef double                                       NT;
typedef leda_rational                                NT;

typedef CGAL::Cartesian<NT>                   K;
typedef CGAL::Partition_traits_2<K>                  Traits;
typedef Traits::Polygon_2                            Polygon;
typedef Traits::Point_2                              Point;
typedef Traits::Segment_2                            Segment;
typedef Polygon::Vertex_const_iterator               Vertex_iterator;

Traits       partition_traits;
int          alg = 0;

void partition(const std::list<Polygon>& input, std::list<Polygon>& output, int tp)
{
  output.clear();  
  std::list<Polygon>::const_iterator it = input.begin();
  
  for(;it != input.end();it++){
     Polygon polygon = *it;
     
     if (! polygon.is_simple()) std::cout << "polygon is not simple !\n";
     else {
       std::cout << "polygon is simple !\n";
     }
  
     switch(tp){
       case 0: { // optimal_convex_partition      
        CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                         polygon.vertices_end(),
                                         std::back_inserter(output),
                                         partition_traits);       
	break;
       }
       case 1: { // greene_approx_convex_partition
        CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(), 
                                         polygon.vertices_end(),
                                         std::back_inserter(output));      
	break;       
       }
       default:{ // approx_convex_partition
        CGAL::approx_convex_partition_2(polygon.vertices_begin(), 
                                         polygon.vertices_end(),
                                         std::back_inserter(output));     
	break;              
       }
     }
  }
}

class geo_poly_partition : public geowin_redraw, public geowin_update<std::list<Polygon>, std::list<Polygon> >
{
public:
  virtual ~geo_poly_partition() {}

  virtual void update(const std::list<Polygon>& L, std::list<Polygon>& output)
  { partition(L,output, alg); }
};

GeoWin* gwin;
geo_scene RES;

void call_back(int choice)
{ alg = choice; RES->update(); gwin->redraw(); }

int main()
{
  geowin_init_default_type((std::list<Polygon>*)0, leda_string("CGALPolygonList"));
 
  std::list<Polygon>  L;

  GeoWin GW("CGAL - polygon partitioning demo");
  gwin = &GW;
  GW.set_bg_color(leda_black);
 
  geo_scene poly_scene= GW.new_scene(L);  
  GW.set_color(poly_scene, leda_red );
  GW.set_line_width(poly_scene, 4);
  GW.set_active_line_width(poly_scene, 4);
  
  geo_poly_partition UPDATE;

  RES = GW.new_scene(UPDATE, poly_scene, leda_string("partition of the polygon")); 
  GW.set_color(RES, leda_green );
  GW.set_line_width(RES, 2);
  GW.set_all_visible(true); 
  
  GW.init_menu();
  leda_list<leda_string> Part;
  Part.append("optimal"); Part.append("greene_approx"); Part.append("approx");
  
  GW.get_window().choice_item("current algorithm:",alg,Part,call_back);
  
  GW.edit(poly_scene);
  
  return 0;  
}

#endif
