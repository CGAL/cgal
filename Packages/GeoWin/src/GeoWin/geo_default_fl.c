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
// file          : src/GeoWin/geo_default_fl.c
// package       : GeoWin (1.2.2)
// revision      : 1.2.2
// revision_date : 30 January 2001 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ============================================================================

#include <LEDA/geowin.h>
#include <LEDA/geowin_init.h>

GEOWIN_BEGIN_NAMESPACE

// post window default input handler ...

void d3_point_window_handler(GeoWin& gw, d3_point& p)
{
   double zv = gw.get_d3_z_value();
   p = d3_point(p.xcoord(),p.ycoord(),zv);
}

// edit functions
void edit_gen_polygon(GeoWin& gw, gen_polygon& gp, int i);
void edit_polygon(GeoWin& gw, polygon& gp, int i);


// insertion of new scenes ...

void insert_new_segment_scene(GeoWin& gw, const list<segment>& LPol)
{
 if (gw.get_algo_create_scenes()) {
  geo_scene sc = GeoScene::new_scene(string("Segments"));
  gw.insert_scene(sc);
  ((GeoEditScene<list<segment> >*)sc)->set_objects(LPol);
  gw.init_and_set_visible(sc);
 } 
}

void insert_new_circle_scene(GeoWin& gw, const list<circle>& LPol)
{
 if (gw.get_algo_create_scenes()) {
  geo_scene sc = GeoScene::new_scene(string("Circles"));
  gw.insert_scene(sc);
  ((GeoEditScene<list<circle> >*)sc)->set_objects(LPol);
  gw.init_and_set_visible(sc);  
 }
}


// get defining points functions ...
void segment_defining_points(const segment& s, list<point>& lp)
{ lp.append(s.source()); lp.append(s.target()); }

void circle_defining_points(const circle& c, list<point>& lp)
{ lp.append(c.point1()); lp.append(c.point2()); lp.append(c.point3()); }

void ray_defining_points(const ray& r, list<point>& lp)
{ lp.append(r.point1()); lp.append(r.point2()); }

void line_defining_points(const line& l, list<point>& lp)
{ lp.append(l.point1()); lp.append(l.point2()); }

void polygon_defining_points(const polygon& p, list<point>& lp)
{ lp=p.vertices(); }

void gen_polygon_defining_points(const gen_polygon& gp, list<point>& lp)
{ lp=gp.vertices(); }

// algorithms ...
void chull_float(GeoWin& gw, list<point>& L);
void triang_float(GeoWin& gw, list<point>& L);
void width_float(GeoWin& gw, list<point>& L);
void dt_float(GeoWin& gw, list<point>& L);
void f_dt_float(GeoWin& gw, list<point>& L);
void mst_float(GeoWin& gw, list<point>& L);
void lec_float(GeoWin& gw, list<point>& L);
void sec_float(GeoWin& gw, list<point>& L);
void all_empty_float(GeoWin& gw, list<point>& L);
void all_encl_float(GeoWin& gw, list<point>& L);
void seg_inter_float(GeoWin& gw, list<segment>& L);
void point_fcn_two_float(GeoWin& gw, list<point>& L);
void point_fcn_three_float(GeoWin& gw, list<point>& L);
void point_fcn_four_float(GeoWin& gw, list<point>& L);
//void check_gen_pol_float(GeoWin& gw, list<gen_polygon>& L);
//void union_gen_pol_float(GeoWin& gw, list<gen_polygon>& L);
//void inter_gen_pol_float(GeoWin& gw, list<gen_polygon>& L);

// input functions/objects ...

void segment_points( list<point>& L, segment s, int n);

void circle_points(list<rat_point>& L, circle C, int n);

class alt_input_d3_point : public GeoInputObject<d3_point>
{
 void operator()(GeoWin& gw, list<d3_point>& LP)
 {
  window& w = gw.get_window();
  point p;
  w >> p;
  
  double x0 = p.xcoord();
  double y0 = p.ycoord();
  double x1, y1;
  w.read_mouse_seg(x0, y0, x1, y1);
  
  LP.append(d3_point(x0,y1,x1-x0));
 }
};

class alt_input_segment : public GeoInputObject<point>
{
int pt_number;
public:
 alt_input_segment() { pt_number= 20; }

 void operator()(GeoWin& gw, list<point>& LP)
 {
  window& w = gw.get_window();
  segment s;
  w >> s;
  segment_points(LP,s,pt_number); 
 }
 
 void options(GeoWin& gw) 
 {
   panel P;
   P.text_item("\\bf\\blue Number of points on segment:");
   P.int_item("number", pt_number, 10, 200);
   P.button(" DONE ", 0);
   gw.open_panel(P);
 }
};

class alt_input_circle : public GeoInputObject<point>
{
int pt_number;
public:
 alt_input_circle() { pt_number= 20; }
 
 void operator()(GeoWin& gw, list<point>& LP)
 {
  window& w = gw.get_window();
  circle c;
  w >> c;
  list<rat_point> LPR;
  circle_points(LPR,c,pt_number); 
  rat_point iter;
  forall(iter,LPR) LP.append(iter.to_float());
 }
 
 void options(GeoWin& gw) 
 {
   panel P;
   P.text_item("\\bf\\blue Number of points on circle:");
   P.int_item("number", pt_number, 10, 200);
   P.button(" DONE ", 0);
   gw.open_panel(P);
 } 
};


class alt_input_square : public GeoInputObject<polygon> {
public:
 void operator()(GeoWin& gw, list<polygon>& LP)
 {
  window& w = gw.get_window();
  gw.message("square : input a segment");
  
  segment RS;
  w >> RS;
  segment RS2 = RS.rotate90(), RS3 = RS2.translate(RS.to_vector());
  point P1 = RS.source(), P2 = RS.target(), P3 = RS2.target(), P4 = RS3.target(); 
  list<point> poly_points;
  poly_points.append(P1);  poly_points.append(P2); poly_points.append(P4); poly_points.append(P3);
  LP.append(polygon(poly_points));
  gw.message("");
 }
};

class alt_input_parallelogram : public GeoInputObject<polygon> {
public:
 void operator()(GeoWin& gw, list<polygon>& LP)
 {
  window& w = gw.get_window();

  segment RS;
  w >> RS; w << RS;
  point P1 = RS.source(), P2 = RS.target();
  double x0=P2.xcoord(),y0=P2.ycoord(), x1, y1;
  w.read_mouse_seg(x0,y0,x1,y1);
  point P3(x1,y1), P4 = P3.translate(-RS.to_vector());
  
  list<point> poly_points;
  poly_points.append(P1);  poly_points.append(P2); poly_points.append(P3); poly_points.append(P4);
  LP.append(polygon(poly_points));
 }
};

/*
class alt_input_gen_poly : public GeoInputObject<gen_polygon>
{
public:
 void operator()(GeoWin& gw, list<gen_polygon>& LG)
 {
  window& w = gw.get_window();
  gw.message("Middle mouse button - finish gen_polygon input");
  
  gen_polygon GP;
  polygon in;

  int but;
  unsigned long t;
  double x,y;
  
  do {
   w.read_event(but, x, y, t);
 
   if (but != MOUSE_BUTTON(2)) {
    w >> in; w << in;
    if (! in.empty() || in.orientation()==0) {
       if (in.orientation()==1) GP = GP.unite(gen_polygon(in));
       else GP = GP.diff(gen_polygon(in.complement()));
    }
   }
  } while (but != MOUSE_BUTTON(2));
  
  gw.message("");
  if (! GP.empty())  LG.append(GP);
 }
};
*/

static alt_input_segment         altseg;
static alt_input_circle          altcirc;
static alt_input_square          altsquare;
static alt_input_parallelogram   altpara;
//static alt_input_gen_poly        altpoly;
static alt_input_d3_point        altd3point;

// ----------------------------------------------
// d3 output functions ...
// ----------------------------------------------

void points_d3(const list<point>& L,d3_window& W, GRAPH<d3_point,int>& H);
void d3_points_d3(const list<d3_point>& L,d3_window& W, GRAPH<d3_point,int>& H);
void segments_d3(const list<segment>& L,d3_window& W, GRAPH<d3_point,int>& H);
void circles_d3(const list<circle>& L,d3_window& W, GRAPH<d3_point,int>& H);
void lines_d3(const list<line>& L,d3_window& W, GRAPH<d3_point,int>& H);
void rays_d3(const list<ray>& L,d3_window& W, GRAPH<d3_point,int>& H);
void poly_d3(const list<polygon>& L,d3_window& W, GRAPH<d3_point,int>& H);
void gen_poly_d3(const list<gen_polygon>& L,d3_window& W, GRAPH<d3_point,int>& H);


window& d3_point_draw(window& w, const d3_point& p,int val)
{
  //cout << val << "\n";
  val = GeoWin::get_projection_mode();
  
  switch(val){
    case 0:  //xy
    { w << point(p.xcoord(),p.ycoord()); break; }
    case 1:  //xz
    { w << point(p.xcoord(),p.zcoord()); break; }    
    case 2:  //yz
    { w << point(p.ycoord(),p.zcoord()); break; }
  }
  return w;
}

#if defined(__KCC)

typedef void (*POINT_D3)(const list<point>&, d3_window&, GRAPH<d3_point,int>&);
typedef void (*SEGMENT_D3)(const list<segment>&, d3_window&, GRAPH<d3_point,int>&);
typedef void (*RAY_D3)(const list<ray>&, d3_window&, GRAPH<d3_point,int>&);
typedef void (*LINE_D3)(const list<line>&, d3_window&, GRAPH<d3_point,int>&);
typedef void (*CIRCLE_D3)(const list<circle>&, d3_window&, GRAPH<d3_point,int>&);
typedef void (*POLYGON_D3)(const list<polygon>&, d3_window&, GRAPH<d3_point,int>&);
typedef void (*GENPOLYGON_D3)(const list<gen_polygon>&, d3_window&, GRAPH<d3_point,int>&);
typedef void (*D3POINT_D3)(const list<d3_point>&, d3_window&, GRAPH<d3_point,int>&);

typedef window& (*D3POINT_DRAW_FCN)(window&, const d3_point&,int);

template class GeoEditScene<list<point> >;
template class GeoBaseScene<list<point> >;
template string geowin_info_fcn(const list<point>&);
template GeoEditScene<list<point> >* make_edit_prototype(list<point>*, string);
template GeoBaseScene<list<point> >* make_base_prototype(list<point>*, string);
template GeoEditScene<list<point> >* geowin_init_default_type(list<point>*, string, POINT_D3); 

template class GeoEditScene<list<segment> >;
template class GeoBaseScene<list<segment> >;
template string geowin_info_fcn(const list<segment>&);
template GeoEditScene<list<segment> >* make_edit_prototype(list<segment>*, string);
template GeoBaseScene<list<segment> >* make_base_prototype(list<segment>*, string);
template GeoEditScene<list<segment> >* geowin_init_default_type(list<segment>*, string, SEGMENT_D3); 

template class GeoEditScene<list<line> >;
template class GeoBaseScene<list<line> >;
template string geowin_info_fcn(const list<line>&);
template GeoEditScene<list<line> >* make_edit_prototype(list<line>*, string);
template GeoBaseScene<list<line> >* make_base_prototype(list<line>*, string);
template GeoEditScene<list<line> >* geowin_init_default_type(list<line>*, string, LINE_D3); 

template class GeoEditScene<list<circle> >;
template class GeoBaseScene<list<circle> >;
template string geowin_info_fcn(const list<circle>&);
template GeoEditScene<list<circle> >* make_edit_prototype(list<circle>*, string);
template GeoBaseScene<list<circle> >* make_base_prototype(list<circle>*, string);
template GeoEditScene<list<circle> >* geowin_init_default_type(list<circle>*, string, CIRCLE_D3); 

template class GeoEditScene<list<polygon> >;
template class GeoBaseScene<list<polygon> >;
template string geowin_info_fcn(const list<polygon>&);
template GeoEditScene<list<polygon> >* make_edit_prototype(list<polygon>*, string);
template GeoBaseScene<list<polygon> >* make_base_prototype(list<polygon>*, string);
template GeoEditScene<list<polygon> >* geowin_init_default_type(list<polygon>*, string, POLYGON_D3); 

template class GeoEditScene<list<gen_polygon> >;
template class GeoBaseScene<list<gen_polygon> >;
template string geowin_info_fcn(const list<gen_polygon>&);
template GeoEditScene<list<gen_polygon> >* make_edit_prototype(list<gen_polygon>*, string);
template GeoBaseScene<list<gen_polygon> >* make_base_prototype(list<gen_polygon>*, string);
template GeoEditScene<list<gen_polygon> >* geowin_init_default_type(list<gen_polygon>*, string, GENPOLYGON_D3); 

template class GeoEditScene<list<d3_point> >;
template class GeoBaseScene<list<d3_point> >;
template string geowin_info_fcn(const list<d3_point>&);
template GeoEditScene<list<d3_point> >* make_edit_prototype(list<d3_point>*, string);
template GeoBaseScene<list<d3_point> >* make_base_prototype(list<d3_point>*, string);
template GeoEditScene<list<d3_point> >* geowin_init_default_type(list<d3_point>*, string, D3POINT_D3);


template void geowin_set_draw_object_fcn(GeoBaseScene<list<d3_point> >*, D3POINT_DRAW_FCN); 


#endif

void init_leda_float_default_types()
{
#if defined INIT_GEOWIN_LEDA_DEFAULT_TYPES

  typedef void (*FU)(...);

  GeoEditScene<list<point> >* esc = geowin_init_default_type((list<point>*)0, string("Points"),points_d3);
  esc->add_alternative_input_object(altseg,  string("Segments"));
  esc->add_alternative_input_object(altcirc, string("Circles"));   
  esc->add_algorithm(chull_float, "Convex hull");
  esc->add_algorithm(triang_float, "Triangulation");
  esc->add_algorithm(width_float, "Width");
  esc->add_algorithm(dt_float, "Delaunay triangulation");
  esc->add_algorithm(f_dt_float, "Furthest point DT");
  esc->add_algorithm(mst_float, "MST");
  esc->add_algorithm(lec_float, "Largest empty circle");
  esc->add_algorithm(sec_float, "Smallest enclosing circle");
  esc->add_algorithm(all_empty_float, "All empty circles");
  esc->add_algorithm(all_encl_float, "All enclosing circles");
  esc->add_algorithm(point_fcn_two_float, "Check (2)");
  esc->add_algorithm(point_fcn_three_float, "Check (3)");
  esc->add_algorithm(point_fcn_four_float, "Check (4)");  
  
  GeoEditScene<list<segment> >* esc2 = geowin_init_default_type((list<segment>*)0,  string("Segments"),segments_d3);
  esc2->set_defining_points_fcn(segment_defining_points);
  esc2->move_point_fcn = (FU) geowin_Translatepoint_seg;
  esc2->add_algorithm(seg_inter_float, "Segment intersection");  

/*
  GeoEditScene<list<ray> >* esc3 = geowin_init_default_type((list<ray>*)0,         string("Rays"), rays_d3);
  esc3->set_defining_points_fcn(ray_defining_points);
  esc3->move_point_fcn = (FU) geowin_Translatepoint_ray;
*/
  
  GeoEditScene<list<line> >* esc4 = geowin_init_default_type((list<line>*)0,        string("Lines"), lines_d3);
  esc4->set_defining_points_fcn(line_defining_points);
  esc4->move_point_fcn = (FU) geowin_Translatepoint_line;
  
  GeoEditScene<list<circle> >* esc5 = geowin_init_default_type((list<circle>*)0,      string("Circles"), circles_d3);
  esc5->set_defining_points_fcn(circle_defining_points);
  esc5->move_point_fcn = (FU) geowin_Translatepoint_circle;

  GeoEditScene<list<polygon> >* esc6 = geowin_init_default_type((list<polygon>*)0, string("SimplePolygons"), poly_d3);
  esc6->add_alternative_input_object(altsquare, string("Square"));
  esc6->add_alternative_input_object(altpara, string("Parallelogram"));
  esc6->set_defining_points_fcn(polygon_defining_points);
  esc6->move_point_fcn = (FU) geowin_Translatepoint_poly;
  esc6->edit_obj_fcn= (FU) edit_polygon;
  
  GeoEditScene<list<gen_polygon> >* esc7 = geowin_init_default_type((list<gen_polygon>*)0, \
  string("GeneralizedPolygons"),gen_poly_d3);
  //esc7->add_alternative_input_object(altpoly, string("Construction")); 
  esc7->set_defining_points_fcn(gen_polygon_defining_points);
  esc7->move_point_fcn = (FU) geowin_Translatepoint_gpoly;  
  //esc7->add_algorithm(check_gen_pol_float, string("Check (1)"));
  //esc7->add_algorithm(union_gen_pol_float, string("Union (2)"));
  //esc7->add_algorithm(inter_gen_pol_float, string("Intersection (2)"));  
  esc7->edit_obj_fcn= (FU) edit_gen_polygon;
  
  
#if defined LEDA_RECTANGLES
  geowin_init_default_type((list<rectangle>*)0,   string("Rectangles"));
#endif

  GeoEditScene<list<d3_point> >* esc8 = geowin_init_default_type((list<d3_point>*)0, \
  string("D3-Points"), d3_points_d3);
  esc8->add_alternative_input_object(altd3point,  string("xyz - Input"));
  // draw fcn ...
  // paras.append("x/y");  paras.append("x/z");  paras.append("y/z");
  geowin_set_draw_object_fcn(esc8, d3_point_draw);
  // esc8->draw_object_parameters = paras;
  esc8->post_window_default_input_handler = d3_point_window_handler;
    
#endif
}

GEOWIN_END_NAMESPACE

