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
// file          : demo/GeoWin/CGAL_delaunay_triang_3.C
//
// ======================================================================

#include <CGAL/basic.h>

#if (!defined(CGAL_USE_LEDA) || (__LEDA__ < 400)) || defined(__KCC)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.0 or higher installed!\n";
 std::cout << "A LEDA version >= 4.0 is required to run GeoWin!\n";
 std::cout << "This demo doensn't work on KCC.\n";
 return 0;
}
#else 

#include <CGAL/Cartesian.h>
#include <CGAL/geowin_support.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef double                                               NT;
typedef CGAL::Cartesian<NT>                                  K;
typedef K::Point_3                                           Point_3;
typedef CGAL::Delaunay_triangulation_3<K>                    Delaunay_3;
typedef std::list<Point_3>                                   Point_3_list;

void show_d3_points(geo_scene sc, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GeoEditScene<Point_3_list>* esc = (GeoEditScene<Point_3_list> *) sc;
 Point_3_list& L= esc->get_objref();

 Point_3 p;
 Point_3_list::const_iterator it = L.begin();

 // input points ...
 for(; it != L.end();++it) { p=*it; H.new_node(convert_to_leda(p)); }

 // compute triangulation...
 Delaunay_3  G_delaunay;
 G_delaunay.clear();
 G_delaunay.insert(L.begin(),L.end());

 // construct graph for visualization ...
 typedef Delaunay_3::Vertex_handle Vertex_handle;
 typedef Delaunay_3::Cell_handle Cell_handle;
 typedef Delaunay_3::Finite_edges_iterator Finite_edges_iterator;

 Vertex_handle v1, v2;
 Cell_handle f;
 int n1, n2;

 Finite_edges_iterator eit = G_delaunay.finite_edges_begin();
 Finite_edges_iterator beyond = G_delaunay.finite_edges_end();

 // edges ...
 for ( ;eit != beyond; ++eit) {
	f  = (*eit).first;
	n1 = (*eit).second;
       	n2 = (*eit).third;
	v1 = f->vertex(n1);
	v2 = f->vertex(n2);
        leda_node n1 = H.new_node(convert_to_leda(v1->point()));
	leda_node n2 = H.new_node(convert_to_leda(v2->point()));
	leda_edge e1 = H.new_edge(n1,n2);
	leda_edge e2 = H.new_edge(n2,n1);
	H.set_reversal(e1,e2);
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
  GeoWin GW("CGAL 3d Delaunay triangulation");
  GW.message("To show the 3d Delaunay triangulation use 'Show d3 output' in Window menu.");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_d3_fcn(my_scene, show_d3_points);

  GW.edit(my_scene);
  
  return 0;  
}

#endif
