// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : src/GeoWin/geo_algo.c
// package       : GeoWin (1.2.2)
// revision      : 1.2.2
// revision_date : 30 January 2001 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ============================================================================


#include <LEDA/geowin.h>
#include <LEDA/rat_geo_alg.h>
#include <LEDA/float_geo_alg.h>

GEOWIN_BEGIN_NAMESPACE

void insert_new_segment_scene(GeoWin& gw, const list<segment>& LPol);
void insert_new_rat_segment_scene(GeoWin& gw, const list<rat_segment>& LPol);
void insert_new_circle_scene(GeoWin& gw, const list<circle>& LPol);
void insert_new_rat_circle_scene(GeoWin& gw, const list<rat_circle>& LPol);

//void insert_new_gen_polygon_scene(GeoWin& gw, const list<gen_polygon>& LPol);
//void insert_new_rat_gen_polygon_scene(GeoWin& gw, const list<rat_gen_polygon>& LPol);


// Convex Hulls 
void chull_float(GeoWin& gw, list<point>& L)
{
/*
  polygon Pl = CONVEX_HULL_POLY(L);
  list<segment> Lseg = Pl.segments();

  window& win = gw.get_window();
  color cold= win.set_color(blue);
  segment siter;
  forall(siter,Lseg) win << siter;  
  win.set_color(cold);
  
  insert_new_segment_scene(gw,Lseg);
*/
}

void chull_rat(GeoWin& gw, list<rat_point>& L)
{
  rat_polygon Pl = CONVEX_HULL_POLY(L);
  list<rat_segment> Lseg = Pl.segments();
    
  window& win = gw.get_window();
  color cold= win.set_color(blue); 
  rat_segment siter;
  forall(siter,Lseg) win << siter;
  win.set_color(cold);

  insert_new_rat_segment_scene(gw,Lseg);  
}

// Triangulations
void triang_float(GeoWin& gw, list<point>& L)
{
  GRAPH<point,int> G;
  TRIANGULATE_POINTS(L,G); 
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  edge e;
  list<segment> LPol;  
  forall_edges(e,G) {
   win << segment(G[source(e)],G[target(e)]);
   LPol.append(segment(G[source(e)],G[target(e)]));
  }
  win.set_color(cold);  
  
  insert_new_segment_scene(gw,LPol);
}

void triang_rat(GeoWin& gw, list<rat_point>& L)
{
  GRAPH<rat_point,int> G;
  TRIANGULATE_POINTS(L,G); 
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  edge e;
  list<rat_segment> LPol;  
  forall_edges(e,G) {
   win << rat_segment(G[source(e)],G[target(e)]);
   LPol.append(rat_segment(G[source(e)],G[target(e)]));
  }
  win.set_color(cold);  

  insert_new_rat_segment_scene(gw,LPol);  
}

// WIDTH
void width_float(GeoWin& gw, list<point>& L)
{
  line l1,l2;
  WIDTH(L,l1,l2);
  window& win = gw.get_window();
  color cold= win.set_color(blue); 
  win << l1 << l2;
  win.set_color(cold);  
}

void width_rat(GeoWin& gw, list<rat_point>& L)
{
  rat_line l1,l2;
  WIDTH(L,l1,l2);
  window& win = gw.get_window();
  color cold= win.set_color(blue); 
  win << l1 << l2;
  win.set_color(cold);  
}
 
// Delaunay triangulation

void dt_float(GeoWin& gw, list<point>& L)
{
  GRAPH<point,int> G;
  DELAUNAY_TRIANG(L,G); 
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  edge e;
  list<segment> LPol;
  forall_edges(e,G) {
   win << segment(G[source(e)],G[target(e)]);   
   LPol.append(segment(G[source(e)],G[target(e)]));
  }
  win.set_color(cold);  
  
  insert_new_segment_scene(gw,LPol);
}

void dt_rat(GeoWin& gw, list<rat_point>& L)
{
  GRAPH<rat_point,int> G;
  DELAUNAY_TRIANG(L,G); 
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  edge e;
  list<rat_segment> LPol;
  forall_edges(e,G) {
   win << rat_segment(G[source(e)],G[target(e)]);
   LPol.append(rat_segment(G[source(e)],G[target(e)]));
  }
  win.set_color(cold);  

  insert_new_rat_segment_scene(gw,LPol);  
}

// Furthest point Delaunay Triangulation

void f_dt_float(GeoWin& gw, list<point>& L)
{
  GRAPH<point,int> G;
  F_DELAUNAY_TRIANG(L,G); 
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  edge e;
  list<segment> LPol;
  forall_edges(e,G) {
   win << segment(G[source(e)],G[target(e)]);
   LPol.append(segment(G[source(e)],G[target(e)]));
  }
  win.set_color(cold);  

  insert_new_segment_scene(gw,LPol);    
}

void f_dt_rat(GeoWin& gw, list<rat_point>& L)
{
  GRAPH<rat_point,int> G;
  F_DELAUNAY_TRIANG(L,G); 
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  edge e;
  list<rat_segment> LPol;
  forall_edges(e,G) {
   win << rat_segment(G[source(e)],G[target(e)]);
   LPol.append(rat_segment(G[source(e)],G[target(e)]));   
  }
  win.set_color(cold);  

  insert_new_rat_segment_scene(gw,LPol);      
}

// MST

void mst_float(GeoWin& gw, list<point>& L)
{
  GRAPH<point,int> G;
  MIN_SPANNING_TREE(L,G); 
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  edge e;
  list<segment> LPol;
  forall_edges(e,G) {
   win << segment(G[source(e)],G[target(e)]);
   LPol.append(segment(G[source(e)],G[target(e)]));
  }
  win.set_color(cold); 
  
  insert_new_segment_scene(gw,LPol);        
}

void mst_rat(GeoWin& gw, list<rat_point>& L)
{
  GRAPH<rat_point,int> G;
  MIN_SPANNING_TREE(L,G); 
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  edge e;
  list<rat_segment> LPol;
  forall_edges(e,G) {
   win << rat_segment(G[source(e)],G[target(e)]);
   LPol.append(rat_segment(G[source(e)],G[target(e)]));
  }
  win.set_color(cold); 
  
  insert_new_rat_segment_scene(gw,LPol);      
}

// LEC

void lec_float(GeoWin& gw, list<point>& L)
{
  circle C = LARGEST_EMPTY_CIRCLE(L);
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  win << C;
  win.set_color(cold);   
  
  list<circle> LPol;
  LPol.append(C);
  insert_new_circle_scene(gw,LPol);    
}

void lec_rat(GeoWin& gw, list<rat_point>& L)
{
  rat_circle C = LARGEST_EMPTY_CIRCLE(L);
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  win << C;
  win.set_color(cold);   
  
  list<rat_circle> LPol;
  LPol.append(C);
  insert_new_rat_circle_scene(gw,LPol);   
}

// SEC

void sec_float(GeoWin& gw, list<point>& L)
{
  circle C = SMALLEST_ENCLOSING_CIRCLE(L);
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  win << C;
  win.set_color(cold);   

  list<circle> LPol;
  LPol.append(C);
  insert_new_circle_scene(gw,LPol);      
}

void sec_rat(GeoWin& gw, list<rat_point>& L)
{
  rat_circle C = SMALLEST_ENCLOSING_CIRCLE(L);
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  win << C;
  win.set_color(cold);   

  list<rat_circle> LPol;
  LPol.append(C);
  insert_new_rat_circle_scene(gw,LPol);   
}

// all empty circles

void all_empty_float(GeoWin& gw, list<point>& L)
{
  list<circle> LC;
  circle C;
  ALL_EMPTY_CIRCLES(L,LC);
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  forall(C,LC) win << C;
  win.set_color(cold);   

  insert_new_circle_scene(gw,LC);   
}

void all_empty_rat(GeoWin& gw, list<rat_point>& L)
{
  list<rat_circle> LC;
  rat_circle C;
  ALL_EMPTY_CIRCLES(L,LC);
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  forall(C,LC) win << C;
  win.set_color(cold);  

  insert_new_rat_circle_scene(gw,LC);   
}

// all enclosing circles

void all_encl_float(GeoWin& gw, list<point>& L)
{
  list<circle> LC;
  circle C;
  ALL_ENCLOSING_CIRCLES((const list<point>&)L,LC);
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  forall(C,LC) win << C;
  win.set_color(cold);   

  insert_new_circle_scene(gw,LC);   
}

void all_encl_rat(GeoWin& gw, list<rat_point>& L)
{
  list<rat_circle> LC;
  rat_circle C;
  ALL_ENCLOSING_CIRCLES((const list<rat_point>&)L,LC);
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  forall(C,LC) win << C;
  win.set_color(cold);   

  insert_new_rat_circle_scene(gw,LC);   
}

// point fcns

void point_fcn_two_float(GeoWin& gw, list<point>& L)
{
 if ( L.size() != 2) return;
 point p1 = L.head();
 point p2 = L[L.succ(L.first())]; 
 string str_p1,str_p2;
 string str_sqrdist, str_xdist, str_ydist;
 string_ostream O1,O2,O3,O4,O5;
 O1 << p1 << ends; str_p1=O1.str();
 O2 << p2 << ends; str_p2=O2.str();
 O3 << p1.sqr_dist(p2) << ends; str_sqrdist=O3.str();
 O4 << p1.xdist(p2) << ends; str_xdist=O4.str();
 O5 << p1.ydist(p2) << ends; str_ydist=O5.str(); 
 panel P;
 P.text_item("point p1:" + str_p1);
 P.text_item("");
 P.text_item("point p2:" + str_p2);
 P.text_item("");
 P.text_item("squared distance:" + str_sqrdist);
 P.text_item("");
 P.text_item("x distance:" + str_xdist);
 P.text_item("");
 P.text_item("y distance:" + str_ydist);
 P.text_item("");
   
 P.button("OK");
 gw.open_panel(P);
}
 
void point_fcn_two_rat(GeoWin& gw, list<rat_point>& L)
{
 if ( L.size() != 2) return;
 rat_point p1 = L.head();
 rat_point p2 = L[L.succ(L.first())];  
 string str_p1,str_p2;
 string str_sqrdist, str_xdist, str_ydist;
 string_ostream O1,O2,O3,O4,O5;
 O1 << p1 << ends; str_p1=O1.str();
 O2 << p2 << ends; str_p2=O2.str();
 O3 << p1.sqr_dist(p2) << ends; str_sqrdist=O3.str();
 O4 << p1.xdist(p2) << ends; str_xdist=O4.str();
 O5 << p1.ydist(p2) << ends; str_ydist=O5.str(); 
 panel P;
 P.text_item("point p1:" + str_p1);
 P.text_item("");
 P.text_item("point p2:" + str_p2);
 P.text_item("");
 P.text_item("squared distance:" + str_sqrdist);
 P.text_item("");
 P.text_item("x distance:" + str_xdist);
 P.text_item("");
 P.text_item("y distance:" + str_ydist);
 P.text_item(""); 
 P.button("OK");
 gw.open_panel(P);
}

void point_fcn_three_float(GeoWin& gw, list<point>& L)
{
 if ( L.size() != 3) return;
 point p1 = L.head();
 point p2 = L[L.succ(L.first())]; 
 point p3 = L[L.get_item(2)];
 string str_p1, str_p2, str_p3;
 string str_ori, str_area;
 string_ostream O1,O2,O3,O4,O5;
 O1 << p1 << ends; str_p1=O1.str();
 O2 << p2 << ends; str_p2=O2.str();
 O3 << p3 << ends; str_p3=O3.str();
 O4 << orientation(p1,p2,p3) << ends; str_ori=O4.str();
 O5 << area(p1,p2,p3) << ends; str_area=O5.str(); 
 panel P;
 P.text_item("point p1:" + str_p1);
 P.text_item("");
 P.text_item("point p2:" + str_p2);
 P.text_item("");
 P.text_item("point p3:" + str_p3);
 P.text_item("");
 P.text_item("orientation(p1,p2,p3):" + str_ori);
 P.text_item("");
 P.text_item("area(p1,p2,p3):" + str_area);
 P.text_item(""); 
 P.button("OK");
 gw.open_panel(P);
}

void point_fcn_three_rat(GeoWin& gw, list<rat_point>& L)
{
 if ( L.size() != 3) return;
 rat_point p1 = L.head();
 rat_point p2 = L[L.succ(L.first())]; 
 rat_point p3 = L[L.get_item(2)];
 string str_p1, str_p2, str_p3;
 string str_ori, str_area;
 string_ostream O1,O2,O3,O4,O5;
 O1 << p1 << ends; str_p1=O1.str();
 O2 << p2 << ends; str_p2=O2.str();
 O3 << p3 << ends; str_p3=O3.str();
 O4 << orientation(p1,p2,p3) << ends; str_ori=O4.str();
 O5 << area(p1,p2,p3) << ends; str_area=O5.str(); 
 panel P;
 P.text_item("point p1:" + str_p1);
 P.text_item("");
 P.text_item("point p2:" + str_p2);
 P.text_item("");
 P.text_item("point p3:" + str_p3);
 P.text_item("");
 P.text_item("orientation(p1,p2,p3):" + str_ori);
 P.text_item("");
 P.text_item("area(p1,p2,p3):" + str_area);
 P.text_item(""); 
 P.button("OK");
 gw.open_panel(P);
}

void point_fcn_four_float(GeoWin& gw, list<point>& L)
{
 if ( L.size() != 4) return;
 point p1 = L.head();
 point p2 = L[L.succ(L.first())]; 
 point p3 = L[L.get_item(2)];
 point p4 = L[L.get_item(3)];
 string str_p1, str_p2, str_p3, str_p4;
 string str_sos;
 string_ostream O1,O2,O3,O4,O5;
 O1 << p1 << ends; str_p1=O1.str();
 O2 << p2 << ends; str_p2=O2.str();
 O3 << p3 << ends; str_p3=O3.str();
 O4 << p4 << ends; str_p4=O4.str();
 O5 << side_of_circle(p1,p2,p3,p4) << ends; str_sos=O5.str(); 
 panel P;
 P.text_item("point p1:" + str_p1);
 P.text_item("");
 P.text_item("point p2:" + str_p2);
 P.text_item("");
 P.text_item("point p3:" + str_p3);
 P.text_item("");
 P.text_item("point p4:" + str_p4);
 P.text_item(""); 
 P.text_item("side of circle(p1,p2,p3,p4):" + str_sos);
 P.text_item("");
 P.button("OK");
 gw.open_panel(P);
}

void point_fcn_four_rat(GeoWin& gw, list<rat_point>& L)
{
 if ( L.size() != 4) return;
 rat_point p1 = L.head();
 rat_point p2 = L[L.succ(L.first())]; 
 rat_point p3 = L[L.get_item(2)];
 rat_point p4 = L[L.get_item(3)];
 string str_p1, str_p2, str_p3, str_p4;
 string str_sos;
 string_ostream O1,O2,O3,O4,O5;
 O1 << p1 << ends; str_p1=O1.str();
 O2 << p2 << ends; str_p2=O2.str();
 O3 << p3 << ends; str_p3=O3.str();
 O4 << p4 << ends; str_p4=O4.str();
 O5 << side_of_circle(p1,p2,p3,p4) << ends; str_sos=O5.str(); 
 panel P;
 P.text_item("point p1:" + str_p1);
 P.text_item("");
 P.text_item("point p2:" + str_p2);
 P.text_item("");
 P.text_item("point p3:" + str_p3);
 P.text_item("");
 P.text_item("point p4:" + str_p4);
 P.text_item(""); 
 P.text_item("side of circle(p1,p2,p3,p4):" + str_sos);
 P.text_item("");
 P.button("OK");
 gw.open_panel(P);
}

// segment intersection

void seg_inter_float(GeoWin& gw, list<segment>& S)
{
  //cout << "seg inter float!\n";
  list<point> LI;
  SEGMENT_INTERSECTION(S,LI);
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  point_style pold = win.set_point_style(disc_point);
  point iter;
  forall(iter,LI) win << iter;
  win.set_color(cold);
  win.set_point_style(pold);
}

void seg_inter_rat(GeoWin& gw, list<rat_segment>& S)
{
  //cout << "seg inter rat!\n";
  list<rat_point> LI;
  SEGMENT_INTERSECTION(S,LI);
  window& win = gw.get_window();
  color cold= win.set_color(blue);
  point_style pold = win.set_point_style(disc_point);
  rat_point iter;
  forall(iter,LI) win << iter;
  win.set_color(cold);
  win.set_point_style(pold);
}


// edit functions ...

// edit functions

// decls.
bool geowin_IntersectsBox(const segment&,  double, double, double, double, bool);
bool geowin_IntersectsBox(const rat_segment&, double, double, double, double, bool);

void edit_gen_polygon(GeoWin& gw, gen_polygon& gp, int i)
{
  list<point> LP = gp.vertices();
  list<segment> LS = gp.edges();

  window& w = gw.get_window();
  gw.message("Left button - delete vertex/insert vertex near edge, right - quit edit");
  double d  = w.pix_to_real(1)*gw.get_mfct();
  int but;
  unsigned long t;
  double x,y;
  int event;
  
  do {
   do event = w.read_event(but, x, y, t);   
   while (event != button_release_event);
   
   double x1=x-d, x2=x+d, y1=y-d, y2=y+d;
   list_item iter,found=NULL;
   int cnt=0,number=0;
   
   if (but == MOUSE_BUTTON(1)){ // left ...
      // delete vertex ...
      forall_items(iter,LP){
        double xw = LP[iter].xcoord(), yw = LP[iter].ycoord();
        if ((xw>=x1) && (xw<=x2) && (yw>=y1) && (yw<=y2)) found=iter;
      }
      if (found != NULL) { 
         LP.del_item(found);
	 gp = gen_polygon(LP);
	 LP = gp.vertices(); LS = gp.edges();
	 gw.get_active_scene()->update_and_redraw();
      }   
      else {
       // insert vertex ...
       forall_items(iter,LS) {
        if (geowin_IntersectsBox(LS[iter],x1,y1,x2,y2,true)) { found=iter; number=cnt; }
	cnt++;
       }
       if (found != NULL) { // insert vertex ...
         list_item it = LP.get_item(number);
	 LP.insert(point(x,y),it);
	 gp = gen_polygon(LP);
	 LP = gp.vertices(); LS = gp.edges();
	 gw.get_active_scene()->update_and_redraw();         
       }
      }
   }

  } while (but != MOUSE_BUTTON(3));
  
  gw.message("");
}


void edit_polygon(GeoWin& gw, polygon& gp, int i)
{
  list<point> LP = gp.vertices();
  list<segment> LS = gp.edges();

  window& w = gw.get_window();
  gw.message("Left button - delete vertex/insert vertex near edge, right - quit edit");
  double d  = w.pix_to_real(1)*gw.get_mfct();
  int but;
  unsigned long t;
  double x,y;
  int event;
  
  do {
   do event = w.read_event(but, x, y, t);   
   while (event != button_release_event);
   
   double x1=x-d, x2=x+d, y1=y-d, y2=y+d;
   list_item iter,found=NULL;
   int cnt=0,number=0;

   if (but == MOUSE_BUTTON(1)){ // left ...
      // delete vertex ...
      forall_items(iter,LP){
        double xw = LP[iter].xcoord(), yw = LP[iter].ycoord();
        if ((xw>=x1) && (xw<=x2) && (yw>=y1) && (yw<=y2)) found=iter;
      }
      if (found != NULL) { 
         LP.del_item(found);
	 gp = polygon(LP);
	 LP = gp.vertices(); LS = gp.edges();
	 gw.get_active_scene()->update_and_redraw();
      }   
      else {
       // insert vertex ...
       forall_items(iter,LS) {
        if (geowin_IntersectsBox(LS[iter],x1,y1,x2,y2,true)) { found=iter; number=cnt; }
	cnt++;
       }
       if (found != NULL) { // insert vertex ...
         list_item it = LP.get_item(number);
	 LP.insert(point(x,y),it);
	 gp = polygon(LP);
	 LP = gp.vertices(); LS = gp.edges();
	 gw.get_active_scene()->update_and_redraw();         
       }
      }
   }   
  } while (but != MOUSE_BUTTON(3));
  
  gw.message("");
}

void edit_rat_gen_polygon(GeoWin& gw, rat_gen_polygon& gp, int i)
{
  list<rat_point> LP = gp.vertices();
  list<rat_segment> LS = gp.edges();

  window& w = gw.get_window();
  gw.message("Left button - delete vertex/insert vertex near edge, right - quit edit");  
  double d  = w.pix_to_real(1)*gw.get_mfct();
  int but;
  unsigned long t;
  double x,y;
  int event;
  
  do {
   do event = w.read_event(but, x, y, t);   
   while (event != button_release_event);

   double x1=x-d, x2=x+d, y1=y-d, y2=y+d;
   list_item iter,found=NULL;
   int cnt=0,number=0;
   
   if (but == MOUSE_BUTTON(1)){ // left ...
      forall_items(iter,LP){
        double xw = LP[iter].to_float().xcoord(), yw = LP[iter].to_float().ycoord();
        if ((xw>=x1) && (xw<=x2) && (yw>=y1) && (yw<=y2)) found=iter;
      }
      if (found != NULL) { // delete a vertex ...
         LP.del_item(found);
	 gp = rat_gen_polygon(LP);
	 LP = gp.vertices(); LS = gp.edges();
	 gw.get_active_scene()->update_and_redraw();
      }
      else {
       forall_items(iter,LS) {
        if (geowin_IntersectsBox(LS[iter],x1,y1,x2,y2,true)) { found=iter; number=cnt; }
	cnt++;
       }
       if (found != NULL) { // insert vertex ...
         list_item it = LP.get_item(number);
	 LP.insert(rat_point(point(x,y)),it);
	 gp = rat_gen_polygon(LP);
	 LP = gp.vertices(); LS = gp.edges();
	 gw.get_active_scene()->update_and_redraw();         
       }
      }
   }
  } while (but != MOUSE_BUTTON(3));
  
  gw.message("");
}


void edit_rat_polygon(GeoWin& gw, rat_polygon& gp, int i)
{
  list<rat_point> LP = gp.vertices();
  list<rat_segment> LS = gp.edges();

  window& w = gw.get_window();
  gw.message("Left button - delete vertex/insert vertex near edge, right - quit edit");
  double d  = w.pix_to_real(1)*gw.get_mfct();
  int but;
  unsigned long t;
  double x,y;
  int event;
  
  do {
   do event = w.read_event(but, x, y, t);   
   while (event != button_release_event);

   double x1=x-d, x2=x+d, y1=y-d, y2=y+d;
   list_item iter,found=NULL;
   int cnt=0,number=0;

   if (but == MOUSE_BUTTON(1)){ // left ...
      forall_items(iter,LP){
        double xw = LP[iter].to_float().xcoord(), yw = LP[iter].to_float().ycoord();
        if ((xw>=x1) && (xw<=x2) && (yw>=y1) && (yw<=y2)) found=iter;
      }
      if (found != NULL) { // delete a vertex ...
         LP.del_item(found);
	 gp = rat_polygon(LP);
	 LP = gp.vertices(); LS = gp.edges();
	 gw.get_active_scene()->update_and_redraw();
      }
      else {
       forall_items(iter,LS) {
        if (geowin_IntersectsBox(LS[iter],x1,y1,x2,y2,true)) { found=iter; number=cnt; }
	cnt++;
       }
       if (found != NULL) { // insert vertex ...
         list_item it = LP.get_item(number);
	 LP.insert(rat_point(point(x,y)),it);
	 gp = rat_polygon(LP);
	 LP = gp.vertices(); LS = gp.edges();
	 gw.get_active_scene()->update_and_redraw();         
       }
      }
   }
  } while (but != MOUSE_BUTTON(3));
  
  gw.message("");
}


// ----------------------------------------------
// d3 output functions ...
// ----------------------------------------------

// float

void points_d3(const list<point>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 point iter;
 forall(iter,L) G.new_node(d3_point(iter.xcoord(),iter.ycoord(),0));
 H.join(G);
}

void d3_points_d3(const list<d3_point>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_point iter;
 forall(iter,L) G.new_node(iter);
 H.join(G);
}

void segments_d3(const list<segment>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 segment iter;
 forall(iter,L) {
   node v1 = G.new_node(d3_point(iter.source().xcoord(),iter.source().ycoord(),0));
   node v2 = G.new_node(d3_point(iter.target().xcoord(),iter.target().ycoord(),0));   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);
 }
 H.join(G);
}

void circle_segments(list<segment>&, circle, int);

void circles_d3(const list<circle>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 list<segment> LS;
 circle C;
 forall(C,L) {
  LS.clear();
  circle_segments(LS,C,30);
  segments_d3(LS,W,H);
 }
}

void lines_d3(const list<line>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 line iter;
 forall(iter,L) {
   point pm((iter.point1().xcoord()+iter.point2().xcoord())/2,(iter.point1().ycoord()+iter.point2().ycoord())/2);
   vector v = iter.point1() - iter.point2();
   v= v * 100;
   point p1=pm+v, p2=pm-v;
   node v1 = G.new_node(d3_point(p1.xcoord(),p1.ycoord(),0));
   node v2 = G.new_node(d3_point(p2.xcoord(),p2.ycoord(),0));   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 }
 H.join(G);
}

void rays_d3(const list<ray>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 ray iter;
 forall(iter,L) {
   vector v = iter.point2() - iter.point1();
   v= v * 100;
   point p1=iter.source(), p2=p1 + v;
   node v1 = G.new_node(d3_point(p1.xcoord(),p1.ycoord(),0));
   node v2 = G.new_node(d3_point(p2.xcoord(),p2.ycoord(),0));   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 } 
 H.join(G);
}


void poly_d3(const list<polygon>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 polygon piter;
 list<segment> LS;
 forall(piter,L) {
   LS = piter.segments();
   segment siter;
   forall(siter,LS) {
    node v1 = G.new_node(d3_point(siter.source().xcoord(),siter.source().ycoord(),0));
    node v2 = G.new_node(d3_point(siter.target().xcoord(),siter.target().ycoord(),0));   
    edge e1 = G.new_edge(v1,v2);
    edge e2 = G.new_edge(v2,v1);
    G.set_reversal(e1,e2);
   }   
 }
 H.join(G); 
}

void gen_poly_d3(const list<gen_polygon>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 gen_polygon piter;
 list<segment> LS;
 forall(piter,L) {
   LS = piter.edges();
   segment siter;
   forall(siter,LS) {
    node v1 = G.new_node(d3_point(siter.source().xcoord(),siter.source().ycoord(),0));
    node v2 = G.new_node(d3_point(siter.target().xcoord(),siter.target().ycoord(),0));   
    edge e1 = G.new_edge(v1,v2);
    edge e2 = G.new_edge(v2,v1);
    G.set_reversal(e1,e2);
   }   
 }
 H.join(G); 
}


// rat

void rat_points_d3(const list<rat_point>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 rat_point iter;
 forall(iter,L) G.new_node(d3_point(iter.to_float().xcoord(),iter.to_float().ycoord(),0));
 H.join(G);
}

void d3_rat_points_d3(const list<d3_rat_point>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_rat_point iter;
 forall(iter,L) G.new_node(iter.to_float());
 H.join(G);
}

void rat_segments_d3(const list<rat_segment>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 rat_segment iter;
 forall(iter,L) {
   node v1 = G.new_node(d3_point(iter.source().to_float().xcoord(),iter.source().to_float().ycoord(),0));
   node v2 = G.new_node(d3_point(iter.target().to_float().xcoord(),iter.target().to_float().ycoord(),0));   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);
 }
 H.join(G);
}

void circle_segments(list<segment>&, circle, int);

void rat_circles_d3(const list<rat_circle>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 list<segment> LS;
 list<rat_segment> LSR;
 rat_circle C;
 forall(C,L) {
  LS.clear(); LSR.clear();
  circle_segments(LS,C.to_float(),30);
  segment iter;
  forall(iter,LS) LSR.append(rat_segment(iter));
  rat_segments_d3(LSR,W,H);
 }
}

void rat_lines_d3(const list<rat_line>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 rat_line r_iter;
 forall(r_iter,L) {
   line iter = r_iter.to_float();
   point pm((iter.point1().xcoord()+iter.point2().xcoord())/2,(iter.point1().ycoord()+iter.point2().ycoord())/2);
   vector v = iter.point1() - iter.point2();
   v= v * 100;
   point p1=pm+v, p2=pm-v;
   node v1 = G.new_node(d3_point(p1.xcoord(),p1.ycoord(),0));
   node v2 = G.new_node(d3_point(p2.xcoord(),p2.ycoord(),0));   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 }
 H.join(G);
}

void rat_rays_d3(const list<rat_ray>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 rat_ray r_iter;
 forall(r_iter,L) {
   ray iter = r_iter.to_float();
   vector v = iter.point2() - iter.point1();
   v= v * 100;
   point p1=iter.source(), p2=p1 + v;
   node v1 = G.new_node(d3_point(p1.xcoord(),p1.ycoord(),0));
   node v2 = G.new_node(d3_point(p2.xcoord(),p2.ycoord(),0));   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 } 
 H.join(G);
}

void rat_poly_d3(const list<rat_polygon>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 rat_polygon piter;
 list<rat_segment> LS;
 forall(piter,L) {
   LS = piter.segments();
   rat_segment siter;
   forall(siter,LS) {
    node v1 = G.new_node(d3_point(siter.source().to_float().xcoord(),siter.source().to_float().ycoord(),0));
    node v2 = G.new_node(d3_point(siter.target().to_float().xcoord(),siter.target().to_float().ycoord(),0));   
    edge e1 = G.new_edge(v1,v2);
    edge e2 = G.new_edge(v2,v1);
    G.set_reversal(e1,e2);
   }   
 }
 H.join(G); 
}

void rat_gen_poly_d3(const list<rat_gen_polygon>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 rat_gen_polygon piter;
 list<rat_segment> LS;
 forall(piter,L) {
   LS = piter.edges();
   rat_segment siter;
   forall(siter,LS) {
    node v1 = G.new_node(d3_point(siter.source().to_float().xcoord(),siter.source().to_float().ycoord(),0));
    node v2 = G.new_node(d3_point(siter.target().to_float().xcoord(),siter.target().to_float().ycoord(),0));   
    edge e1 = G.new_edge(v1,v2);
    edge e2 = G.new_edge(v2,v1);
    G.set_reversal(e1,e2);
   }   
 }
 H.join(G); 
}




GEOWIN_END_NAMESPACE










