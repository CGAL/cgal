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
// file          : demo/GeoWin/CGAL_bops_polyholes_2.C
//
// ======================================================================

// we show how to use the CGAL Nef polyhedra to represent polygons with holes 
// and how to perform boolean operations on them

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

#include <LEDA/basic.h>

#if defined(LEDA_NAMESPACE)
#else
#ifndef LEDA_NAMESPACE_NAME
#define LEDA_NAMESPACE_NAME
#endif
#endif

#include <CGAL/leda_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/Nef_2/PM_visualizor.h>
#include <CGAL/IO/Nef_polyhedron_2_Window_stream.h>
#include <CGAL/geowin_support.h>

typedef CGAL::Extended_homogeneous<leda_integer>              EK;
typedef CGAL::Nef_polyhedron_2<EK>                            Nef;
typedef Nef::Point                                            EPoint;

// types for Nef polyhedron exploration
typedef Nef::Explorer                                         Explorer;
typedef Explorer::Vertex_const_handle                         Vertex_const_handle;
typedef Explorer::Halfedge_const_handle                       Halfedge_const_handle;
typedef Explorer::Halfedge_around_face_const_circulator       Halfedge_around_face_const_circulator;
typedef Explorer::Face_const_handle                           Face_const_handle;
typedef Explorer::Vertex_const_iterator                       Vertex_const_iterator;
typedef Explorer::Halfedge_const_iterator                     Halfedge_const_iterator;
typedef Explorer::Face_const_iterator                         Face_const_iterator;

// the polygons for Nef construction ....
typedef CGAL::Cartesian<double>                               K;
typedef K::Point_2                                            Point;
typedef K::Segment_2                                          Segment;
typedef CGAL::Polygon_traits_2<K>                             PTraits;
typedef CGAL::Polygon_2<PTraits,std::list<Point> >            Polygon;

// drawing routine using the visualizor
// see Nef Windowstream ...

#define USING(t) typedef typename PMCDEC::t t

namespace CGAL {

template <typename PMCDEC>
class PM_MyBooleanColor 
{
  USING(Vertex_const_handle);   
  USING(Halfedge_const_handle); 
  USING(Face_const_handle);
  USING(Mark);
  
  Color vertex_cl, hedge_cl, face_cl;
  int my_width;
  
public:

  PM_MyBooleanColor(Color vc, Color hc, Color fc, int w) : vertex_cl(vc), hedge_cl(hc), face_cl(fc), my_width(w)
  { }

  Color color(Vertex_const_handle, const Mark& m) const
  { return ( m ? vertex_cl : CGAL::Color(190,190,190) ); }
  
  int width(Vertex_const_handle, const Mark& m) const
  { return 1; }
  
  Color color(Halfedge_const_handle, const Mark& m) const
  { return face_cl; } 
  
  int width(Halfedge_const_handle, const Mark& m) const
  { return my_width; }
  
  Color color(Face_const_handle, const Mark& m) const
  { return ( m ? face_cl : CGAL::WHITE ); }
};


template <typename T>
void draw_nef(CGAL::Window_stream& ws, const Nef_polyhedron_2<T>& P, Color cl, int wd)
{
  typedef Nef_polyhedron_2<T>                       Polyhedron;
  typedef typename T::RT                            RT;
  typedef typename T::Standard_RT                   Standard_RT;
  typedef typename Polyhedron::Topological_explorer TExplorer;

  typedef CGAL::PM_MyBooleanColor<TExplorer>        MyColor;
  typedef CGAL::PM_visualizor<TExplorer,T,MyColor>  Visualizor;

  Explorer D = P.explorer();
  const T& E = Nef_polyhedron_2<T>::EK;

  Standard_RT frame_radius = frame_default;
  E.determine_frame_radius(D.points_begin(),D.points_end(),frame_radius);
  RT::set_R(frame_radius);
  
  MyColor colors(CGAL::BLACK, CGAL::BLACK, cl, wd);
  T kernel;
  
  Visualizor PMV(ws,D,kernel,colors); 
   
  // draw segments underlying halfedges: 
  Halfedge_const_iterator hit = D.halfedges_begin(), hend = D.halfedges_end(); 
  for (; hit != hend; hit++) { if (! D.is_frame_edge(hit)) PMV.draw(hit); }

  // draw points underlying vertices:
  Vertex_const_iterator vit, vend = D.vertices_end();
  for (vit = D.vertices_begin(); vit != vend; ++vit) PMV.draw(vit);
}

} // end of CGAL namespace

Nef create_nef_polyhedron(Polygon& act)
{
  std::list<EPoint> pts;
  Polygon::Vertex_const_iterator it= act.vertices_begin();
  Polygon::Vertex_const_iterator st= act.vertices_end();  
  
  for(;it != st; it++){
    double xc = (*it).x();
    double yc = (*it).y();
    leda_rat_point rp(leda_point(xc,yc));
    pts.push_back(EPoint(rp.X(),rp.Y(),rp.W()));
  } 
  return Nef(pts.begin(), pts.end());
}

// what kind of boolean operation ???
int op_kind = 0;

class geo_nef_construction : 
            public LEDA_NAMESPACE_NAME::geowin_update<std::list<Polygon>, std::list<Segment> >,
            public LEDA_NAMESPACE_NAME::geowin_redraw
{
public:
 // this is the constructed polygon with holes
 Nef construction_result; 
 bool empty;
 CGAL::Color inner_color;
 
 geo_nef_construction(CGAL::Color interior) : empty(true), inner_color(interior) { }

 void update(const std::list<Polygon>& L, std::list<Segment>& Sl)
 {
  empty = true;
  construction_result = Nef();
  if (L.size() < 1) return;
  
  std::list<Polygon>::const_iterator poly_iter = L.begin();
  
  Polygon act = *poly_iter; 
  construction_result = create_nef_polyhedron(act);
  
  // cut the holes into the resulting polygon ...
  poly_iter++;
  for(;poly_iter!=L.end();poly_iter++){
     act = *poly_iter; 
     construction_result = construction_result.difference(create_nef_polyhedron(act));  
  }
  empty = false;
 }
 
 bool is_empty() const { return empty; }
 Nef  get_nef()  const { return construction_result; }
 
 void draw(leda_window& W,leda_color c1,leda_color c2,double x1,double y1,double x2,double y2)
 { if (! empty) CGAL::draw_nef(W, construction_result,inner_color,1); }
};

// Boolean Operation ....

class geo_nef : public LEDA_NAMESPACE_NAME::geowin_update<std::list<Polygon>, std::list<Segment> >,
                public LEDA_NAMESPACE_NAME::geowin_redraw
{
public:
 Nef bop_result; 
 bool empty;
 geo_nef_construction* C1;
 geo_nef_construction* C2;
 
 geo_nef(geo_nef_construction* constr1,geo_nef_construction* constr2) : empty(true), C1(constr1), C2(constr2) 
 { }

 // we don't need any containers; instead we use the 2 constructed nefs
 void update(const std::list<Polygon>&, std::list<Segment>&)
 {
  std::cout << "update ...\n";
  empty= true;
  if ( (!C1->is_empty()) && (!C2->is_empty())){
    // perform the boolean operation ...
    Nef N1 = C1->get_nef();
    Nef N2 = C2->get_nef();
    
    switch (op_kind) {
     case 0:  { std::cout << "intersection...\n"; bop_result = N1.intersection(N2);  break; }
     case 1:  { std::cout << "union...\n"; bop_result = N1.join(N2);  break; } 
     case 2:  { std::cout << "difference...\n"; bop_result = N1.difference(N2);  break; }           
    }
    empty = false;
  }
 }
 
 void draw(leda_window& W,leda_color c1,leda_color c2,double x1,double y1,double x2,double y2)
 { if (! empty) CGAL::draw_nef(W, bop_result, CGAL::RED,4); }
};



LEDA_NAMESPACE_NAME::GeoWin* gwin;
LEDA_NAMESPACE_NAME::geo_scene result;
LEDA_NAMESPACE_NAME::geo_scene holes1, holes2;
  

void call_back(int choice)
{ op_kind = choice; result->update(); gwin->redraw(); }

int main()
{
  LEDA_NAMESPACE_NAME::geowin_init_default_type((std::list<Polygon>*)0, leda_string("CGALPolygonList"));
   
  std::list<Polygon> L1;
  std::list<Polygon> L2;

  LEDA_NAMESPACE_NAME::GeoWin GW("Boolean operations on polygons with holes");
  gwin = &GW;
#if  __LEDA__ < 430
  GW.add_help_text(leda_string("Polyholes_2"));
#else
  GW.add_special_help_text(leda_string("Polyholes_2"),true);
#endif  
  
  // we use these two scenes to construct polygons with holes ...
  LEDA_NAMESPACE_NAME::geo_scene my_scene1= GW.new_scene(L1);  
  LEDA_NAMESPACE_NAME::geo_scene my_scene2= GW.new_scene(L2);  
  
  // these two scenes hold the resulting polygons ...
  geo_nef_construction C1(CGAL::GREEN), C2(CGAL::YELLOW); 
  holes1  = GW.new_scene(C1,C1,my_scene1,leda_string("Polygon with holes - 1")); 
  holes2  = GW.new_scene(C2,C2,my_scene2,leda_string("Polygon with holes - 2"));   
  
  geo_nef update_obj(&C1,&C2); 
  result  = GW.new_scene(update_obj,update_obj,my_scene1,leda_string("Result of boolean operation")); 
  GW.add_dependence(my_scene2, result);
  
  leda_list<leda_string> ops;
  ops.append("intersection");
  ops.append("union");
  ops.append("difference");
  
  GW.init_menu();
  GW.get_window().choice_item("Operation:",op_kind,ops,call_back);
 
  GW.set_all_visible(true); 

  GW.edit(my_scene1);
  
  return 0;  
}

#endif
