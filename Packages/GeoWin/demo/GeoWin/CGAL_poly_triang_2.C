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
// file          : demo/GeoWin/CGAL_poly_triang_2.C
//
// ======================================================================

// triangulating polygons using CGAL ...

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
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/geowin_support.h>
#include <CGAL/leda_rational.h>
#include <map>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

#if !defined(_MSC_VER)
typedef leda_rational coord_type;
#else
typedef double coord_type;
#endif

typedef CGAL::Cartesian<coord_type>                             K;
typedef K::Point_2                                              Point;
typedef K::Segment_2                                            Segment;
typedef CGAL::Constrained_Delaunay_triangulation_2<K>           CT;

typedef CT::Constraint                                          Constraint;
typedef CT::Edge                                                Edge;
typedef CT::Vertex_handle                                       Vertex_handle;
typedef CT::Edge_iterator                                       Edge_iterator;
typedef CT::Edge_circulator                                     Edge_circulator;
typedef CT::Vertex_circulator                                   Vertex_circulator;

//polygon typedefs
typedef CGAL::Polygon_traits_2<K>                               PTraits;
typedef CGAL::Polygon_2<PTraits,std::list<Point> >              Polygon;


//mark edges outside ...
void handle_inc_edges(CT& ct,
                      std::list<Vertex_handle>& vhl,
		      std::map<Edge,int>& out_edges)
{
  std::list<Vertex_handle>::const_iterator vhit = vhl.begin();
  for(;vhit != vhl.end();vhit++){
  
    std::list<Vertex_handle>::const_iterator vhelp;
    if (vhit == vhl.begin()) vhelp = (--vhl.end());
    else { vhelp = vhit; vhelp--; }
  
    Vertex_handle vpred= *vhelp;
    Vertex_handle vact = *vhit;
    
    vhelp = vhit; vhelp++;
    if (vhelp == vhl.end()) vhelp = vhl.begin();
    
    Vertex_handle vsucc= *vhelp;
    
    Edge_circulator ec = ct.incident_edges(vact);
    
    bool pred_found = false;
    bool succ_found = false;
    
    do {
      if (! ct.is_infinite(ec)) {
         Edge ea = *ec;
         Vertex_handle vstart = ea.first->vertex(ct.ccw(ea.second));
	 Vertex_handle vend   = ea.first->vertex(ct.cw(ea.second));	 
	 
	 if (vstart==vpred) pred_found = true;
	 if (vstart==vsucc && pred_found) succ_found=true;
	 
	 if (pred_found) {
	   if (! (vstart==vpred || vstart==vsucc)) out_edges[ea] = 1;
	 }
      }
      ec++;
    }
    while (! succ_found);
  }
}


class poly_triang : public geowin_update<std::list<Polygon>, std::list<Segment> >
{
public:
 void update(const std::list<Polygon>& Li, std::list<Segment>& Sl)
 {
  Sl.clear();       
  
  // iterate on the polygons ...
  std::list<Polygon>::const_iterator pit = Li.begin();
  
  for(;pit != Li.end(); pit++) { 
    Polygon pact = *pit;

    // polygon must be simple ...
    if ((! pact.is_simple()) || (pact.size() < 3)) continue;
    
    //if it is not ccw, reverse it ...
    CGAL::Orientation ori = pact.orientation();
    
    if (ori==CGAL::CLOCKWISE) pact.reverse_orientation();
    
    CT ct;
    Polygon::Edge_const_iterator peit = pact.edges_begin();
    std::list<Vertex_handle>  vhl;
    
    for(;peit != pact.edges_end();peit++){
       Point src = (*peit).source();
       Point trg = (*peit).target();
       
       Vertex_handle vh1 = ct.insert(src);
       Vertex_handle vh2 = ct.insert(trg);
       
       // now insert the constraint ...
       ct.insert(vh1, vh2);
 
       vhl.push_back(vh1);
    }
    // now mark all edges outside the polygon ...
    std::map<Edge,int> outside;
    
    handle_inc_edges(ct, vhl, outside);
     
    Edge_iterator eit = ct.edges_begin();
    Edge_iterator beyond = ct.edges_end();   
    Edge eact;

    while (eit != beyond) {
       eact = *eit;     
       if (outside.count(eact) == 0){ 
         Sl.push_back(ct.segment(eact)); 
       }              
       ++eit;  
    } 
  }
 }
};

int main()
{
  geowin_init_default_type((std::list<Polygon>*)0, leda_string("CGALPolygonList"));

  std::list<Polygon> PL;
  GeoWin GW("CGAL polygon triangulation demo");
  GW.message("Please input simple polygons. The demo triangulates them.");
 
  geo_scene my_scene= GW.new_scene(PL); 
  GW.set_color(my_scene, leda_blue); 
  GW.set_fill_color(my_scene, leda_invisible);
  GW.set_active_line_width(my_scene, 3);

  poly_triang PT;
  geo_scene res2 = GW.new_scene(PT, my_scene, leda_string("Polygon triangulation"));
  GW.set_color(res2, leda_red);
  GW.set_visible(res2, true);
  
  GW.edit(my_scene);
  
  return 0;  
}

#endif
