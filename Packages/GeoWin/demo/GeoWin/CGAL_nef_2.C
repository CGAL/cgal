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
// file          : demo/GeoWin/CGAL_nef_2.C
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
public:

  PM_MyBooleanColor(Color vc, Color hc, Color fc) : vertex_cl(vc), hedge_cl(hc), face_cl(fc)
  { }

  Color color(Vertex_const_handle, const Mark& m) const
  { return ( m ? vertex_cl : CGAL::Color(190,190,190) ); }
  
  int width(Vertex_const_handle, const Mark& m) const
  { return 3; }
  
  Color color(Halfedge_const_handle, const Mark& m) const
  { return ( m ? hedge_cl : CGAL::Color(190,190,190) ); }
  
  int width(Halfedge_const_handle, const Mark& m) const
  { return 3; }
  
  Color color(Face_const_handle, const Mark& m) const
  { return ( m ? face_cl : CGAL::WHITE ); }
};


template <typename T>
void draw_nef(CGAL::Window_stream& ws, const Nef_polyhedron_2<T>& P)
{
  typedef Nef_polyhedron_2<T>                       Polyhedron;
  typedef typename T::RT                            RT;
  typedef typename T::Standard_RT                   Standard_RT;
  typedef typename Polyhedron::Topological_explorer TExplorer;

  typedef CGAL::PM_MyBooleanColor<TExplorer>        MyColor;
  typedef CGAL::PM_visualizor<TExplorer,T,MyColor>  Visualizor;

  TExplorer D = P.explorer();
  const T& E = Nef_polyhedron_2<T>::EK;

  Standard_RT frame_radius = frame_default;
  E.determine_frame_radius(D.points_begin(),D.points_end(),frame_radius);
  RT::set_R(frame_radius);
  
  MyColor colors(CGAL::RED, CGAL::RED, CGAL::RED);
  T kernel;
  
  Visualizor PMV(ws,D,kernel,colors); 
  PMV.draw_map();
}

}

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

class geo_nef : public LEDA_NAMESPACE_NAME::geowin_update<std::list<Polygon>, std::list<Segment> >,
                public LEDA_NAMESPACE_NAME::geowin_redraw
{
public:
 Nef bop_result; 
 Polygon poly1, poly2;
 bool empty;
 
 geo_nef() : empty(true) { }

 void update(const std::list<Polygon>& L, std::list<Segment>& Sl)
 {
  Sl.clear();
  empty = true;

  bop_result = Nef();
  poly1 = Polygon(); poly2 = Polygon();
  
  if (L.size() < 2) return;
  
  // build two Nef polyhedra ...
  std::list<Polygon>::const_iterator poly_iter = L.begin();  
  Polygon act = *poly_iter; poly1 = act;
  Nef N1 = create_nef_polyhedron(act);
  poly_iter++;  
  act = *poly_iter; poly2 = act;
  Nef N2 = create_nef_polyhedron(act); 
  
  //perform bop
  switch (op_kind) {
     case 0:  { std::cout << "intersection...\n"; bop_result = N1.intersection(N2);  break; }
     case 1:  { std::cout << "union...\n"; bop_result = N1.join(N2);  break; } 
     case 2:  { std::cout << "difference...\n"; bop_result = N1.difference(N2);  break; }
     case 3:  { std::cout << "symmetric difference...\n"; bop_result = N1.symmetric_difference(N2);  break; }           
  }
  empty = false;
 }
 
 void draw(leda_window& W,leda_color c1,leda_color c2,double x1,double y1,double x2,double y2)
 {  
   if (! empty) CGAL::draw_nef(W, bop_result);   
   W.set_color(c1);
   W << poly1; W << poly2;
 }
};

LEDA_NAMESPACE_NAME::GeoWin* gwin;
LEDA_NAMESPACE_NAME::geo_scene result;  

void call_back(int choice)
{ op_kind = choice; result->update(); gwin->redraw(); }

int main()
{
  LEDA_NAMESPACE_NAME::geowin_init_default_type((std::list<Polygon>*)0, leda_string("CGALPolygonList"));
   
  std::list<Polygon> L;

  LEDA_NAMESPACE_NAME::GeoWin GW("Boolean operations on 2d nef polyhedra");
  gwin = &GW;
#if  __LEDA__ < 430
  GW.add_help_text(leda_string("Nef_2"));
#else
  GW.add_special_help_text(leda_string("Nef_2"),true);
#endif

  geo_nef update_obj;
  
  LEDA_NAMESPACE_NAME::geo_scene my_scene= GW.new_scene(L);  
  
  result  = GW.new_scene(update_obj,update_obj,my_scene,leda_string("Boolean operations")); 
  GW.set_visible(result,true);
  
  leda_list<leda_string> ops;
  ops.append("intersection");
  ops.append("union");
  ops.append("difference");
  ops.append("symmetric difference"); 
  
  GW.message("Input two simple polygons!"); 
  GW.init_menu();
  GW.get_window().choice_item("Operation:",op_kind,ops,call_back);
 
  GW.edit(my_scene);
  
  return 0;  
}

#endif
