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
// file          : demo/GeoWin/CGAL_convex_hull_3.C
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
#include <CGAL/leda_rational.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/predicates_on_points_3.h>


typedef CGAL::Cartesian<leda_rational>            R;
typedef CGAL::Convex_hull_traits_3<R>             Traits;
typedef Traits::Polyhedron_3                      Polyhedron_3;
typedef Polyhedron_3::Vertex_handle               Vertex_handle;
typedef Polyhedron_3::Halfedge_iterator           Halfedge_iterator;


CGAL::Point_3<CGAL::Cartesian<leda_rational> > conversion_to_rat(CGAL::Point_3<CGAL::Cartesian<double> > pin)
{
 double x = pin.x();
 double y = pin.y();
 double z = pin.z();
 leda_d3_point pf(x,y,z);
 leda_d3_rat_point prat(leda_integer(x*100000),leda_integer(y*100000),leda_integer(z*100000),leda_integer(100000));
 return CGAL::Point_3<CGAL::Cartesian<leda_rational> >(prat.xcoord(),prat.ycoord(),prat.zcoord());
}

CGAL::Point_3<CGAL::Cartesian<double> > conversion_to_float(CGAL::Point_3<CGAL::Cartesian<leda_rational> > pin)
{
 leda_rational x = pin.x();
 leda_rational y = pin.y();
 leda_rational z = pin.z();
 leda_d3_rat_point prat(x,y,z);
 leda_d3_point pf = prat.to_d3_point();
 return CGAL::Point_3<CGAL::Cartesian<double> >(pf.xcoord(),pf.ycoord(),pf.zcoord());
}


static void show_d3_points(geo_scene sc, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GeoEditScene< CGALPoint_3_list >* esc = (GeoEditScene< CGALPoint_3_list > *) sc;
 CGALPoint_3_list& L= esc->get_objref();

 CGALPoint_3 p;
 CGALPoint_3_list::const_iterator it = L.begin();

 for(; it != L.end();++it) { p=*it; H.new_node(convert_to_leda(p)); }

 std::list<CGAL::Point_3<CGAL::Cartesian<leda_rational> > > L2;
 // convert list with FP values to rational values ...
 it = L.begin();
 for(; it != L.end();++it) {
    CGAL::Point_3<CGAL::Cartesian<leda_rational> > rpt = conversion_to_rat(*it);
    L2.push_back(rpt);
 }
  
 CGAL::Object ch_object;

 // compute convex hull 
 CGAL::convex_hull_3(L2.begin(), L2.end(), ch_object);

 // visualize polyhedron...
 Polyhedron_3           Pol;
 CGAL::Segment_3<R>     segment;
 CGAL::Triangle_3<R>    triangle;
 CGAL::Point_3<R>       point;
 
 if (assign(Pol, ch_object)) {
  Halfedge_iterator hit = Pol.halfedges_begin();
  for (; hit != Pol.halfedges_end(); hit++) {
    Vertex_handle v1 = hit->vertex();
    Vertex_handle v2 = hit->opposite()->vertex();
    CGALPoint_3 ps = conversion_to_float(v1->point());
    CGALPoint_3 pt = conversion_to_float(v2->point());

    leda_node n1= H.new_node(convert_to_leda(ps));
    leda_node n2= H.new_node(convert_to_leda(pt));
    leda_edge e1= H.new_edge(n1,n2), e2= H.new_edge(n2,n1);
    H.set_reversal(e1,e2);
  }
 }
 
 else if ( assign(segment, ch_object) ){
    CGALPoint_3 ps = conversion_to_float(segment.source());
    CGALPoint_3 pt = conversion_to_float(segment.target()); 
 
    leda_node n1= H.new_node(convert_to_leda(ps));
    leda_node n2= H.new_node(convert_to_leda(pt));
    leda_edge e1= H.new_edge(n1,n2), e2= H.new_edge(n2,n1);
    H.set_reversal(e1,e2);    
 }
 else if (assign(triangle, ch_object) ){
    CGALPoint_3 ps = conversion_to_float(triangle.vertex(0));
    CGALPoint_3 pt = conversion_to_float(triangle.vertex(1));
    CGALPoint_3 pu = conversion_to_float(triangle.vertex(2)); 
    
    leda_node n1= H.new_node(convert_to_leda(ps));
    leda_node n2= H.new_node(convert_to_leda(pt));
    leda_node n3= H.new_node(convert_to_leda(pu));
    leda_edge e1= H.new_edge(n1,n2), e2= H.new_edge(n2,n1);
    H.set_reversal(e1,e2);  
    leda_edge e3= H.new_edge(n2,n3), e4= H.new_edge(n3,n2);
    H.set_reversal(e3,e4);  
    leda_edge e5= H.new_edge(n3,n1), e6= H.new_edge(n1,n3);
    H.set_reversal(e5,e6);              
 }
 else if (assign(point, ch_object) ){
    CGALPoint_3 p = conversion_to_float(point); 
    H.new_node(convert_to_leda(p));
 } 
 
 
 leda_node_array<leda_vector> pos(H);
 leda_node v;
 forall_nodes(v,H) pos[v] = H[v].to_vector();
 W.init(pos); 
}

int main()
{
  geowin_init_default_type((CGALPoint_3_list*)0, leda_string("CGALPoint_3_List"));
 
  CGALPoint_3_list L;
  GeoWin GW("Convex hull in 3d");
  GW.add_help_text(leda_string("CGAL_convex_hull_3"));
  GW.message("To show the convex hull use 'Show d3 output' in Window menu.");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_d3_fcn(my_scene, show_d3_points);

  GW.edit(my_scene);
  
  return 0;  
}

#endif
