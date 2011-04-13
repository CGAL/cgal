// ======================================================================
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
// release_date  : 2001, January 30
//
// file          : src/GeoWin/LEDA/geowin_init_d3.h
// package       : GeoWin (1.2.2)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.2.2
// revision_date : 30 January 2000 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ======================================================================



#ifndef LEDA_GEOWIN_INIT_D3_H
#define LEDA_GEOWIN_INIT_D3_H

#if !defined(LEDA_ROOT_INCL_ID)
#define LEDA_ROOT_INCL_ID 400955
#include <LEDA/REDEFINE_NAMES.h>
#endif

#include <LEDA/geowin.h>
#include <LEDA/geowin_init.h>
#include <LEDA/d3_float_kernel.h>
#include <LEDA/d3_rat_kernel.h>
#include <LEDA/d3_simplex.h>
#include <LEDA/d3_rat_simplex.h>
#include <LEDA/d3_ray.h>
#include <LEDA/d3_rat_ray.h>
#include <LEDA/d3_hull.h>

#if (__LEDA__ >= 420)
#include<LEDA/d3_triangle.h>
#include<LEDA/d3_rat_triangle.h>
#endif


//generator routine to sphere d3 output
void generate_points_on_sphere(const d3_sphere& Sph, list<d3_point>& pts)
{
  pts.clear();
  d3_point ct = Sph.center();
  double r  = Sph.radius();
  double ymin = Sph.center().ycoord() - 0.95*r;
  double ymax = Sph.center().ycoord() + 0.95*r;
  double yakt = ymin, ystep = 0.12*r;
  double xmin,xmax,xakt;
  circle C(point(ct.xcoord(),ct.ycoord()),r);
  circle C2;
  double xstep;
  //cout << C << "\n";
  
  while(yakt <= ymax){
    segment S(point(ct.xcoord()-r-100.0,yakt), point(ct.xcoord()+r+100.0,yakt));
    list<point> res = C.intersection(S);
    point p1 = res.pop();
    point p2 = res.pop();
    if (p1.xcoord() < p2.xcoord()) { xmin = p1.xcoord(); xmax =p2.xcoord(); }
    else { xmin = p2.xcoord(); xmax =p1.xcoord(); }
    xstep=(xmax-xmin)*0.13;    
    pts.append(d3_point(xmin,yakt,ct.zcoord())); pts.append(d3_point(xmax,yakt,ct.zcoord()));        
    xakt=xmin+xstep;
    C2 = circle(point(ct.xcoord(),ct.zcoord()),ct.xcoord()-xmin); 
    
    while (xakt <= xmax){
      // get z -value
      segment S2(point(xakt,ct.zcoord()-r-100.0), point(xakt,ct.zcoord()+r+100.0));
      list<point> res2 = C2.intersection(S2);
      if (res2.size()>1) {
       point p3 = res2.pop();
       point p4 = res2.pop();
      
       pts.append(d3_point(xakt,yakt,p3.ycoord())); 
       pts.append(d3_point(xakt,yakt,p4.ycoord()));  
      }
      xakt = xakt + xstep;
    }
    yakt = yakt + ystep;
  }
  
  pts.append(d3_point(ct.xcoord(),ct.ycoord()-r,ct.zcoord()));
  pts.append(d3_point(ct.xcoord(),ct.ycoord()+r,ct.zcoord()));
}

bool geowin_IntersectsBox(const point& obj, double x1,double y1,double x2, double y2,bool f);
bool geowin_IntersectsBox(const circle& obj, double x1,double y1,double x2, double y2,bool f);
void geowin_BoundingBox(const circle& obj, double& x1, double& x2,double& y1, double& y2);

d3_rat_point leda_convert_to(const d3_point& p)
{
  d3_rat_point pr;
  double x=p.xcoord(), y=p.ycoord(), z=p.zcoord();
  pr = d3_rat_point(integer(x*100000),integer(y*100000), integer(z*100000),100000);  
  return pr;
}

d3_rat_line leda_convert_to(const d3_line& l)
{
  d3_rat_line lr(leda_convert_to(l.point1()),leda_convert_to(l.point2()));
  return lr;
}

d3_rat_segment leda_convert_to(const d3_segment& s)
{
  d3_rat_segment sr(leda_convert_to(s.source()),leda_convert_to(s.target()));
  return sr;
}

d3_rat_ray leda_convert_to(const d3_ray& r)
{
  d3_rat_ray rr(leda_convert_to(r.point1()),leda_convert_to(r.point2()));
  return rr;
}

d3_rat_simplex leda_convert_to(const d3_simplex& s)
{
  d3_rat_simplex sr(leda_convert_to(s.point1()),leda_convert_to(s.point2()), \
                    leda_convert_to(s.point3()),leda_convert_to(s.point4()));
  return sr;
}

d3_rat_triangle leda_convert_to(const d3_triangle& s)
{
  d3_rat_triangle sr(leda_convert_to(s.point1()),leda_convert_to(s.point2()), \
                    leda_convert_to(s.point3()));
  return sr;
}

/*
#if !defined(__SUNPRO_CC) 
d3_rat_polygon leda_convert_to(const d3_polygon& s)
{
  list<d3_rat_point> LPR;
  list<d3_point> LP = s.vertices();
  d3_point pi;
  forall(pi, LP) LPR.append(leda_convert_to(pi));

  d3_rat_polygon sr(LPR);
  return sr;
}
#endif
*/

d3_rat_sphere leda_convert_to(const d3_sphere& s)
{
  d3_rat_sphere sr(leda_convert_to(s.point1()),leda_convert_to(s.point2()), \
                    leda_convert_to(s.point3()),leda_convert_to(s.point4()));
  return sr;
}


// d3_line / d3_rat_line

void d3_lines_d3(const list<d3_line>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_line iter;
 forall(iter,L) {
   d3_point pm((iter.point1().xcoord()+iter.point2().xcoord())/2, \
               (iter.point1().ycoord()+iter.point2().ycoord())/2, \
	       (iter.point1().zcoord()+iter.point2().zcoord())/2);
   vector v = iter.point1() - iter.point2();
   v= v * 100;
   d3_point p1=pm+v, p2=pm-v;
   node v1 = G.new_node(p1);
   node v2 = G.new_node(p2);   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 }
 H.join(G);
}

void d3_rat_lines_d3(const list<d3_rat_line>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_rat_line iter;
 forall(iter,L) {
   d3_rat_point pm((iter.point1().xcoord()+iter.point2().xcoord())/2, \
               (iter.point1().ycoord()+iter.point2().ycoord())/2, \
	       (iter.point1().zcoord()+iter.point2().zcoord())/2);
   rat_vector v = iter.point1() - iter.point2();
   v= integer(100) * v ;
   d3_rat_point p1=pm+v, p2=pm-v;
   node v1 = G.new_node(p1.to_float());
   node v2 = G.new_node(p2.to_float());   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 }
 H.join(G);
}

ps_file& operator<<(ps_file& F,const d3_line& l)
{
   line lh;
   bool b = l.project_xy(lh);
   if (b) F << lh;
   else F << point(l.point1().xcoord(), l.point1().ycoord());
   
   return F;
}

ps_file& operator<<(ps_file& F,const d3_rat_line& l)
{
   rat_line lh;
   bool b = l.project_xy(lh);
   if (b) F << lh;
   else F << rat_point(l.point1().xcoord(), l.point1().ycoord());
   
   return F;
}

window& operator << (window& w, const d3_line& l)
{
   line lh;
   bool b = l.project_xy(lh);
   if (b) w << lh;
   else w << point(l.point1().xcoord(), l.point1().ycoord());
   return w;
}

window& operator << (window& w, const d3_rat_line& l)
{
   rat_line lh;
   bool b = l.project_xy(lh);
   if (b) w << lh;
   else w << rat_point(l.point1().xcoord(), l.point1().ycoord());
   return w;
}

window& d3_line_draw(window& w, const d3_line& l,int val)
{
  val = GeoWin::get_projection_mode();
  bool b;
  line lh;
  
  switch(val){
    case 0:  //xy
    { b = l.project_xy(lh);
      if (b) w << lh;
      else w << point(l.point1().xcoord(), l.point1().ycoord());
      break; }
    case 1:  //xz
    { b = l.project_xz(lh);
      if (b) w << lh;
      else w << point(l.point1().xcoord(), l.point1().zcoord());
      break; }    
    default:  //yz
    { b = l.project_yz(lh);
      if (b) w << lh;
      else w << point(l.point1().ycoord(), l.point1().zcoord()); 
      break; }
  }
  return w;
}

window& d3_rat_line_draw(window& w, const d3_rat_line& l,int val)
{
  return d3_line_draw(w,l.to_float(),val);
}


window& operator >> (window& w, d3_line& obj)
{
   line l;
   w >> l;
   d3_point p1(l.point1().xcoord(), l.point1().ycoord(),0);
   d3_point p2(l.point2().xcoord(), l.point2().ycoord(),0);
   obj = d3_line(p1,p2);
   return w;
}

window& operator >> (window& w, d3_rat_line& obj)
{
   rat_line l;
   w >> l;
   d3_rat_point p1(l.point1().xcoord(), l.point1().ycoord(), 0);
   d3_rat_point p2(l.point2().xcoord(), l.point2().ycoord(), 0);
   obj = d3_rat_line(p1,p2);
   return w;
}

bool geowin_IntersectsBox(const d3_line& obj, double x1,double y1,double x2, double y2,bool f)
{
  int val = GeoWin::get_projection_mode();
  
  segment s1(x1, y1, x2, y1);
  segment s2(x1, y2, x2, y2);
  segment s3(x1, y1, x1, y2);
  segment s4(x2, y1, x2, y2);
  point p;
  
  line lh;
  bool b;
  
  switch(val){
   case 0:  { b = obj.project_xy(lh); break; }
   case 1:  { b = obj.project_xz(lh); break; }
   default:{ b = obj.project_yz(lh); break; }
  }
  
  if (b) {
    return 
    (  lh.intersection(s1,p) || lh.intersection(s2,p) ||
       lh.intersection(s3,p) || lh.intersection(s4,p)  );
  }
  else {
    point p;
    switch(val){
     case 0: {  p = point(obj.point1().xcoord(),obj.point1().ycoord()); break; }
     case 1: {  p = point(obj.point1().xcoord(),obj.point1().zcoord()); break; }
     default: {  p = point(obj.point1().ycoord(),obj.point1().zcoord()); break; }
    }
    return ((p.xcoord() < x2) && (p.xcoord() > x1) && (p.ycoord() < y2) && (p.ycoord() > y1));
  }
}

void geowin_BoundingBox(const d3_line& obj, double& x1, double& x2,double& y1, double& y2)
{  
  int val = GeoWin::get_projection_mode();
  point p1,p2;
  
  switch(val){
   case 0: {
    p1=point(obj.point1().xcoord(),obj.point1().ycoord());
    p2=point(obj.point2().xcoord(),obj.point2().ycoord());
    break;
   }
   case 1: {
    p1=point(obj.point1().xcoord(),obj.point1().zcoord());
    p2=point(obj.point2().xcoord(),obj.point2().zcoord());
    break;
   }
   default: {
    p1=point(obj.point1().ycoord(),obj.point1().zcoord());
    p2=point(obj.point2().ycoord(),obj.point2().zcoord());
    break;
   }   
  }
  
  segment s(p1,p2);

  if (s.xcoord1()<s.xcoord2()) { x1=s.xcoord1(); x2=s.xcoord2(); }
  else { x2=s.xcoord1(); x1=s.xcoord2(); }

  if (s.ycoord1()<s.ycoord2()) { y1=s.ycoord1(); y2=s.ycoord2(); }
  else { y2=s.ycoord1(); y1=s.ycoord2(); }
}

void geowin_Translate(d3_line& obj, double dx, double dy)
{
  int val = GeoWin::get_projection_mode();
  d3_line obnew;
  
  switch(val){
   case 0:   { obnew = obj.translate(dx,dy,0); break; }
   case 1:   { obnew = obj.translate(dx,0,dy); break; }
   default: { obnew = obj.translate(0,dx,dy); break; }
  }
  
  obj = obnew;
}

void geowin_Rotate(d3_line& obj, double dx, double dy,double a){ }

bool geowin_IntersectsBox(const d3_rat_line& obj, double x1,double y1,double x2, double y2,bool f)
{ return geowin_IntersectsBox(obj.to_float(),x1,y1,x2,y2,f); }

void geowin_BoundingBox(const d3_rat_line& obj, double& x1, double& x2,double& y1, double& y2)
{ geowin_BoundingBox(obj.to_float(),x1,x2,y1,y2); }

void geowin_Translate(d3_rat_line& obj, double dx, double dy)
{ 
 int val = GeoWin::get_projection_mode();
 switch(val){
  case 0: { obj = leda_convert_to(obj.to_float().translate(dx, dy, 0)); break; }
  case 1: { obj = leda_convert_to(obj.to_float().translate(dx, 0, dy)); break; }
  default: { obj = leda_convert_to(obj.to_float().translate(0,dx,dy)); break; }  
 }
}

void geowin_Rotate(d3_rat_line& obj, double dx, double dy,double a){ }


void geowin_generate_objects(GeoWin& gw, list<d3_line>& L){ }

void geowin_generate_objects(GeoWin& gw, list<d3_rat_line>& L){ }



// d3_segment / d3_rat_segment
// d3 functions

void d3_segments_d3(const list<d3_segment>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_segment iter;
 forall(iter,L) {
   node v1 = G.new_node(iter.source());
   node v2 = G.new_node(iter.target());   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);
 }
 H.join(G);
}

void d3_rat_segments_d3(const list<d3_rat_segment>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_rat_segment iter;
 forall(iter,L) {
   node v1 = G.new_node(iter.source().to_float());
   node v2 = G.new_node(iter.target().to_float());   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);
 }
 H.join(G);
}


ps_file& operator<<(ps_file& F,const d3_segment& s)
{
   segment sh = s.project_xy();
   F << sh; 
   return F;
}

ps_file& operator<<(ps_file& F,const d3_rat_segment& s)
{
   rat_segment sh = s.project_xy();
   F << sh;
   return F;
}

window& operator << (window& w, const d3_segment& obj)
{
   w << obj.project_xy();
   return w;
}

window& operator << (window& w, const d3_rat_segment& obj)
{
   w << obj.project_xy();
   return w;
}

window& d3_segment_draw(window& w, const d3_segment& obj,int val)
{
  val = GeoWin::get_projection_mode();
  
  switch(val){
    case 0:  //xy
    { w << obj.project_xy(); break; }
    case 1:  //xz
    { w << obj.project_xz(); break; }    
    default:  //yz
    { w << obj.project_yz(); break; }
  }
  return w;
}

window& d3_rat_segment_draw(window& w, const d3_rat_segment& obj,int val)
{
  return d3_segment_draw(w,obj.to_float(),val);
}

window& operator >> (window& w, d3_segment& obj)
{
   segment s;
   w >> s;
   d3_point p1(s.start().xcoord(), s.start().ycoord(), 0);
   d3_point p2(s.end().xcoord(), s.end().ycoord(), 0);
   obj = d3_segment(p1,p2);
   return w;
}

window& operator >> (window& w, d3_rat_segment& obj)
{
   rat_segment s;
   w >> s;
   d3_rat_point p1(s.start().xcoord(), s.start().ycoord(), 0);
   d3_rat_point p2(s.end().xcoord(), s.end().ycoord(), 0);
   obj = d3_rat_segment(p1,p2);
   return w;
}

bool geowin_IntersectsBox(const d3_segment& obj, double x1,double y1,double x2, double y2,bool f)
{
  int val = GeoWin::get_projection_mode();
  segment s1(x1, y1, x2, y1);
  segment s2(x1, y2, x2, y2);
  segment s3(x1, y1, x1, y2);
  segment s4(x2, y1, x2, y2);
  point p;
    
  segment s;
  
  switch(val){
   case 0: { s = obj.project_xy(); break; }
   case 1: { s = obj.project_xz(); break; }
   default: { s = obj.project_yz(); break; }
  }

  return 
    ( geowin_IntersectsBox(s.start(), x1, y1, x2, y2,f) ||
      geowin_IntersectsBox(s.end(), x1, y1, x2, y2,f)   ||
      s1.intersection(s,p) || s2.intersection(s,p) ||
      s3.intersection(s,p) || s4.intersection(s,p) );   
}

void geowin_BoundingBox(const d3_segment& obj, double& x1, double& x2,double& y1, double& y2)
{  
  int val = GeoWin::get_projection_mode();
  segment s;
  
  switch(val){
   case 0: { s = obj.project_xy(); break; }
   case 1: { s = obj.project_xz(); break; }
   default: { s = obj.project_yz(); break; }
  } 
  
  if (s.xcoord1()<s.xcoord2()) { x1=s.xcoord1(); x2=s.xcoord2(); }
  else { x2=s.xcoord1(); x1=s.xcoord2(); }

  if (s.ycoord1()<s.ycoord2()) { y1=s.ycoord1(); y2=s.ycoord2(); }
  else { y2=s.ycoord1(); y1=s.ycoord2(); }
}

void geowin_Translate(d3_segment& obj, double dx, double dy)
{
  int val = GeoWin::get_projection_mode();
  d3_segment obnew;
  
  switch(val){
   case 0:   { obnew = obj.translate(dx,dy,0); break; }
   case 1:   { obnew = obj.translate(dx,0,dy); break; }
   default: { obnew = obj.translate(0,dx,dy); break; }
  } 
  obj = obnew;
}

void geowin_Rotate(d3_segment& obj, double dx, double dy,double a){ }

bool geowin_IntersectsBox(const d3_rat_segment& obj, double x1,double y1,double x2, double y2,bool f)
{ return geowin_IntersectsBox(obj.to_float(),x1,y1,x2,y2,f); }

void geowin_BoundingBox(const d3_rat_segment& obj, double& x1, double& x2,double& y1, double& y2)
{ geowin_BoundingBox(obj.to_float(),x1,x2,y1,y2); }

void geowin_Translate(d3_rat_segment& obj, double dx, double dy)
{ 
 int val = GeoWin::get_projection_mode();
 switch(val){
  case 0: { obj = leda_convert_to(obj.to_float().translate(dx, dy, 0)); break; }
  case 1: { obj = leda_convert_to(obj.to_float().translate(dx, 0, dy)); break; }
  default: { obj = leda_convert_to(obj.to_float().translate(0,dx,dy)); break; }  
 }
}


void geowin_Rotate(d3_rat_segment& obj, double dx, double dy,double a){ }

void geowin_generate_objects(GeoWin& gw, list<d3_segment>& L){ }

void geowin_generate_objects(GeoWin& gw, list<d3_rat_segment>& L){ }



// d3_ray / d3_rat_ray

void d3_rays_d3(const list<d3_ray>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_ray iter;
 forall(iter,L) {
   vector v = iter.point2() - iter.point1();
   v= v * 100;
   d3_point p1=iter.source(), p2=p1 + v;
   node v1 = G.new_node(p1);
   node v2 = G.new_node(p2);   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 } 
 H.join(G);
}

void d3_rat_rays_d3(const list<d3_rat_ray>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_rat_ray iter;
 forall(iter,L) {
   rat_vector v = iter.point2() - iter.point1();
   v= integer(100) * v;
   d3_rat_point p1=iter.source(), p2=p1 + v;
   node v1 = G.new_node(p1.to_float());
   node v2 = G.new_node(p2.to_float());   
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 } 
 H.join(G);
}

ps_file& operator<<(ps_file& F,const d3_ray& r)
{
   ray rh;
   bool b = r.project_xy(rh);
   if (b) F << rh;
   else F << point(r.point1().xcoord(), r.point1().ycoord()); 
   return F;
}

ps_file& operator<<(ps_file& F,const d3_rat_ray& r)
{
   rat_ray rh;
   bool b = r.project_xy(rh);
   if (b) F << rh;
   else F << rat_point(r.point1().xcoord(), r.point1().ycoord());
   
   return F;
}

window& operator << (window& w, const d3_ray& r)
{
   ray rh;
   bool b = r.project_xy(rh);
   if (b) w << rh;
   else w << point(r.point1().xcoord(), r.point1().ycoord());
   return w;
}

window& operator << (window& w, const d3_rat_ray& r)
{
   rat_ray rh;
   bool b = r.project_xy(rh);
   if (b) w << rh;
   else w << rat_point(r.point1().xcoord(), r.point1().ycoord());
   return w;
}

window& d3_ray_draw(window& w, const d3_ray& obj,int val)
{
  val = GeoWin::get_projection_mode();
  bool b;
  ray lh;
  
  switch(val){
    case 0:  //xy
    { b = obj.project_xy(lh);
      if (b) w << lh;
      else w << point(obj.point1().xcoord(), obj.point1().ycoord());
      break; }
    case 1:  //xz
    { b = obj.project_xz(lh);
      if (b) w << lh;
      else w << point(obj.point1().xcoord(), obj.point1().zcoord());
      break; }    
    default:  //yz
    { b = obj.project_yz(lh);
      if (b) w << lh;
      else w << point(obj.point1().ycoord(), obj.point1().zcoord()); 
      break; }
  }
  return w;
}

window& d3_rat_ray_draw(window& w, const d3_rat_ray& obj,int val)
{
  return d3_ray_draw(w,obj.to_float(),val);
}


window& operator >> (window& w, d3_ray& obj)
{
   ray r;
   w >> r;
   d3_point p1(r.point1().xcoord(), r.point1().ycoord(), 0);
   d3_point p2(r.point2().xcoord(), r.point2().ycoord(), 0);
   obj = d3_ray(p1,p2);
   return w;
}

window& operator >> (window& w, d3_rat_ray& obj)
{
   rat_ray r;
   w >> r;
   d3_rat_point p1(r.point1().xcoord(), r.point1().ycoord(), 0);
   d3_rat_point p2(r.point2().xcoord(), r.point2().ycoord(), 0);
   obj = d3_rat_ray(p1,p2);
   return w;
}

bool geowin_IntersectsBox(const d3_ray& obj, double x1,double y1,double x2, double y2,bool f)
{
  int val = GeoWin::get_projection_mode();
  segment s1(x1, y1, x2, y1);
  segment s2(x1, y2, x2, y2);
  segment s3(x1, y1, x1, y2);
  segment s4(x2, y1, x2, y2);
  point p;
  
  ray rh;
  bool b;

  switch(val){
   case 0:  { b = obj.project_xy(rh); break; }
   case 1:  { b = obj.project_xz(rh); break; }
   default:{ b = obj.project_yz(rh); break; }
  } 
  
  if (b) {
    return 
    (  rh.intersection(s1,p) || rh.intersection(s2,p) ||
       rh.intersection(s3,p) || rh.intersection(s4,p)  );
  }
  else {
    point p;
    switch(val){
     case 0: {  p = point(obj.point1().xcoord(),obj.point1().ycoord()); break; }
     case 1: {  p = point(obj.point1().xcoord(),obj.point1().zcoord()); break; }
     default: {  p = point(obj.point1().ycoord(),obj.point1().zcoord()); break; }
    } 
    return ((p.xcoord() < x2) && (p.xcoord() > x1) && (p.ycoord() < y2) && (p.ycoord() > y1));
  }
}

void geowin_BoundingBox(const d3_ray& obj, double& x1, double& x2,double& y1, double& y2)
{  
  int val = GeoWin::get_projection_mode();
  point p1,p2;
  
  switch(val){
   case 0: {
    p1=point(obj.point1().xcoord(),obj.point1().ycoord());
    p2=point(obj.point2().xcoord(),obj.point2().ycoord());
    break;
   }
   case 1: {
    p1=point(obj.point1().xcoord(),obj.point1().zcoord());
    p2=point(obj.point2().xcoord(),obj.point2().zcoord());
    break;
   }
   default: {
    p1=point(obj.point1().ycoord(),obj.point1().zcoord());
    p2=point(obj.point2().ycoord(),obj.point2().zcoord());
    break;
   }   
  }
  segment s(p1,p2);

  if (s.xcoord1()<s.xcoord2()) { x1=s.xcoord1(); x2=s.xcoord2(); }
  else { x2=s.xcoord1(); x1=s.xcoord2(); }

  if (s.ycoord1()<s.ycoord2()) { y1=s.ycoord1(); y2=s.ycoord2(); }
  else { y2=s.ycoord1(); y1=s.ycoord2(); }
}

void geowin_Translate(d3_ray& obj, double dx, double dy)
{
  d3_ray obnew; 
  int val = GeoWin::get_projection_mode();
  
  switch(val){
   case 0:   { obnew = obj.translate(dx,dy,0); break; }
   case 1:   { obnew = obj.translate(dx,0,dy); break; }
   default: { obnew = obj.translate(0,dx,dy); break; }
  }  
  obj = obnew;
}

void geowin_Rotate(d3_ray& obj, double dx, double dy,double a){ }

bool geowin_IntersectsBox(const d3_rat_ray& obj, double x1,double y1,double x2, double y2,bool f)
{ return geowin_IntersectsBox(obj.to_float(),x1,y1,x2,y2,f); }

void geowin_BoundingBox(const d3_rat_ray& obj, double& x1, double& x2,double& y1, double& y2)
{ geowin_BoundingBox(obj.to_float(),x1,x2,y1,y2); }

void geowin_Translate(d3_rat_ray& obj, double dx, double dy)
{
 int val = GeoWin::get_projection_mode();
 switch(val){
  case 0: { obj = leda_convert_to(obj.to_float().translate(dx, dy, 0)); break; }
  case 1: { obj = leda_convert_to(obj.to_float().translate(dx, 0, dy)); break; }
  default: { obj = leda_convert_to(obj.to_float().translate(0,dx,dy)); break; }  
 }
}

void geowin_Rotate(d3_rat_ray& obj, double dx, double dy,double a){ }


void geowin_generate_objects(GeoWin& gw, list<d3_ray>& L){ }

void geowin_generate_objects(GeoWin& gw, list<d3_rat_ray>& L){ }


// d3_simplex / d3_rat_simplex

// d3_simplex / d3_rat_simplex
// d3 functions

void d3_simplex_d3(const list<d3_simplex>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_simplex iter;
 forall(iter,L) {
   node v1 = G.new_node(iter.point1());
   node v2 = G.new_node(iter.point2());  
   node v3 = G.new_node(iter.point3());
   node v4 = G.new_node(iter.point4()); 
      
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v1,v4);
   edge e3 = G.new_edge(v1,v3);
   
   edge e4 = G.new_edge(v2,v3);
   edge e5 = G.new_edge(v2,v4);
   edge e6 = G.new_edge(v2,v1);  
   
   edge e7 = G.new_edge(v3,v1);
   edge e8 = G.new_edge(v3,v4);
   edge e9 = G.new_edge(v3,v2);  

   edge e10 = G.new_edge(v4,v2);
   edge e11 = G.new_edge(v4,v3);
   edge e12 = G.new_edge(v4,v1);
   
   // set reversals ...
   G.set_reversal(e1,e6);
   G.set_reversal(e2,e12);
   G.set_reversal(e3,e7);
   G.set_reversal(e4,e9);
   G.set_reversal(e5,e10);
   G.set_reversal(e8,e11);
 }
 H.join(G);
}

void d3_rat_simplex_d3(const list<d3_rat_simplex>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_rat_simplex iter;
 forall(iter,L) {
   node v1 = G.new_node(iter.point1().to_float());
   node v2 = G.new_node(iter.point2().to_float());  
   node v3 = G.new_node(iter.point3().to_float());
   node v4 = G.new_node(iter.point4().to_float()); 
      
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v1,v4);
   edge e3 = G.new_edge(v1,v3);
   
   edge e4 = G.new_edge(v2,v3);
   edge e5 = G.new_edge(v2,v4);
   edge e6 = G.new_edge(v2,v1);  
   
   edge e7 = G.new_edge(v3,v1);
   edge e8 = G.new_edge(v3,v4);
   edge e9 = G.new_edge(v3,v2);  

   edge e10 = G.new_edge(v4,v2);
   edge e11 = G.new_edge(v4,v3);
   edge e12 = G.new_edge(v4,v1);
   
   // set reversals ...
   G.set_reversal(e1,e6);
   G.set_reversal(e2,e12);
   G.set_reversal(e3,e7);
   G.set_reversal(e4,e9);
   G.set_reversal(e5,e10);
   G.set_reversal(e8,e11); 
 }
 H.join(G);
}

ps_file& operator<<(ps_file& F,const d3_simplex& r)
{
   point p1(r.point1().xcoord(), r.point1().ycoord());
   point p2(r.point2().xcoord(), r.point2().ycoord());
   point p3(r.point3().xcoord(), r.point3().ycoord());
   point p4(r.point4().xcoord(), r.point4().ycoord());

   F << segment(p1,p2); F << segment(p1,p3); F << segment(p1,p4);
   F << segment(p2,p3); F << segment(p2,p4); F << segment(p3,p4);
      
   return F;
}

ps_file& operator<<(ps_file& F,const d3_rat_simplex& r)
{
   point p1(r.point1().to_float().xcoord(), r.point1().to_float().ycoord());
   point p2(r.point2().to_float().xcoord(), r.point2().to_float().ycoord());
   point p3(r.point3().to_float().xcoord(), r.point3().to_float().ycoord());
   point p4(r.point4().to_float().xcoord(), r.point4().to_float().ycoord());

   F << segment(p1,p2); F << segment(p1,p3); F << segment(p1,p4);
   F << segment(p2,p3); F << segment(p2,p4); F << segment(p3,p4);
      
   return F;
}

window& operator << (window& w, const d3_simplex& r)
{
   point p1(r.point1().xcoord(), r.point1().ycoord());
   point p2(r.point2().xcoord(), r.point2().ycoord());
   point p3(r.point3().xcoord(), r.point3().ycoord());
   point p4(r.point4().xcoord(), r.point4().ycoord());

   w << segment(p1,p2); w << segment(p1,p3); w << segment(p1,p4);
   w << segment(p2,p3); w << segment(p2,p4); w << segment(p3,p4);
      
   return w;
}

window& operator << (window& w, const d3_rat_simplex& r)
{
   point p1(r.point1().to_float().xcoord(), r.point1().to_float().ycoord());
   point p2(r.point2().to_float().xcoord(), r.point2().to_float().ycoord());
   point p3(r.point3().to_float().xcoord(), r.point3().to_float().ycoord());
   point p4(r.point4().to_float().xcoord(), r.point4().to_float().ycoord());

   w << segment(p1,p2); w << segment(p1,p3); w << segment(p1,p4);
   w << segment(p2,p3); w << segment(p2,p4); w << segment(p3,p4);
      
   return w;
}

window& d3_simplex_draw(window& w, const d3_simplex& obj,int val)
{
  val = GeoWin::get_projection_mode();
  point p1,p2,p3,p4;

  switch(val){
   case 0:{
    p1=point(obj.point1().xcoord(), obj.point1().ycoord());
    p2=point(obj.point2().xcoord(), obj.point2().ycoord());
    p3=point(obj.point3().xcoord(), obj.point3().ycoord());
    p4=point(obj.point4().xcoord(), obj.point4().ycoord()); break;
   }   
   case 1:{
    p1=point(obj.point1().xcoord(), obj.point1().zcoord());
    p2=point(obj.point2().xcoord(), obj.point2().zcoord());
    p3=point(obj.point3().xcoord(), obj.point3().zcoord());
    p4=point(obj.point4().xcoord(), obj.point4().zcoord()); break;
   }   
   default:{
    p1=point(obj.point1().ycoord(), obj.point1().zcoord());
    p2=point(obj.point2().ycoord(), obj.point2().zcoord());
    p3=point(obj.point3().ycoord(), obj.point3().zcoord());
    p4=point(obj.point4().ycoord(), obj.point4().zcoord()); break;
   }   
  }
  
  w << segment(p1,p2); w << segment(p1,p3); w << segment(p1,p4);
  w << segment(p2,p3); w << segment(p2,p4); w << segment(p3,p4);
  
  return w;
}

window& d3_rat_simplex_draw(window& w, const d3_rat_simplex& obj,int val)
{
  return d3_simplex_draw(w,obj.to_float(),val);
}

window& operator >> (window& w, d3_simplex& obj)
{
   d3_point p1,p2,p3,p4;
   point p;
   w >> p;
   double x,y;
   
   p1 = d3_point(p.xcoord(),p.ycoord(),0);
   w.read_mouse_seg(p.xcoord(),p.ycoord(),x,y);
   w << segment(p.xcoord(),p.ycoord(),x,y);
   p2 = d3_point(x,y,0);
   w.read_mouse_seg(p2.xcoord(),p2.ycoord(),x,y);
   p3 = d3_point(x,y,0);
   w << segment(p3.xcoord(),p3.ycoord(),p1.xcoord(),p1.ycoord());
   w << segment(p3.xcoord(),p3.ycoord(),p2.xcoord(),p2.ycoord());   
   w << segment(p2.xcoord(),p2.ycoord(),p1.xcoord(),p1.ycoord());
   w.read_mouse_seg(p3.xcoord(),p3.ycoord(),x,y);
   p4 = d3_point(x,y,0);
   
   obj = d3_simplex(p1,p2,p3,p4);
   return w;
}

window& operator >> (window& w, d3_rat_simplex& obj)
{
   d3_simplex s;
   w >> s;
   obj = leda_convert_to(s);
   return w;
}

bool geowin_IntersectsBox(const d3_simplex& obj, double x1,double y1,double x2, double y2,bool f)
{
  int val = GeoWin::get_projection_mode();
  point p1,p2,p3,p4;

  switch(val){
   case 0:{
    p1=point(obj.point1().xcoord(), obj.point1().ycoord());
    p2=point(obj.point2().xcoord(), obj.point2().ycoord());
    p3=point(obj.point3().xcoord(), obj.point3().ycoord());
    p4=point(obj.point4().xcoord(), obj.point4().ycoord()); break;
   }   
   case 1:{
    p1=point(obj.point1().xcoord(), obj.point1().zcoord());
    p2=point(obj.point2().xcoord(), obj.point2().zcoord());
    p3=point(obj.point3().xcoord(), obj.point3().zcoord());
    p4=point(obj.point4().xcoord(), obj.point4().zcoord()); break;
   }   
   default:{
    p1=point(obj.point1().ycoord(), obj.point1().zcoord());
    p2=point(obj.point2().ycoord(), obj.point2().zcoord());
    p3=point(obj.point3().ycoord(), obj.point3().zcoord());
    p4=point(obj.point4().ycoord(), obj.point4().zcoord()); break;
   }   
  }
    
  rectangle R(x1,y1,x2,y2);
  list<point> IL;
  IL = R.intersection(segment(p1,p2));
  if (! IL.empty()) return true;
  IL = R.intersection(segment(p1,p3));
  if (! IL.empty()) return true;  
  IL = R.intersection(segment(p1,p4));
  if (! IL.empty()) return true;  
  IL = R.intersection(segment(p2,p3));
  if (! IL.empty()) return true;
  IL = R.intersection(segment(p2,p4));
  if (! IL.empty()) return true;  
  IL = R.intersection(segment(p3,p4));
  if (! IL.empty()) return true;  
  return false;
}

void geowin_BoundingBox(const d3_simplex& obj, double& x1, double& x2,double& y1, double& y2)
{ 
  int val = GeoWin::get_projection_mode();
  point p1,p2,p3,p4;
  switch(val){
   case 0:{
    p1=point(obj.point1().xcoord(), obj.point1().ycoord());
    p2=point(obj.point2().xcoord(), obj.point2().ycoord());
    p3=point(obj.point3().xcoord(), obj.point3().ycoord());
    p4=point(obj.point4().xcoord(), obj.point4().ycoord()); break;
   }   
   case 1:{
    p1=point(obj.point1().xcoord(), obj.point1().zcoord());
    p2=point(obj.point2().xcoord(), obj.point2().zcoord());
    p3=point(obj.point3().xcoord(), obj.point3().zcoord());
    p4=point(obj.point4().xcoord(), obj.point4().zcoord()); break;
   }   
   default:{
    p1=point(obj.point1().ycoord(), obj.point1().zcoord());
    p2=point(obj.point2().ycoord(), obj.point2().zcoord());
    p3=point(obj.point3().ycoord(), obj.point3().zcoord());
    p4=point(obj.point4().ycoord(), obj.point4().zcoord()); break;
   }   
  } 
  
  x1=p1.xcoord(); x2=p1.xcoord(); y1=p1.ycoord(); y2=p1.ycoord();
  // p2
  if (p2.xcoord()<x1) x1=p2.xcoord();
  if (p2.xcoord()>x2) x2=p2.xcoord(); 
  if (p2.ycoord()<y1) y1=p2.ycoord();
  if (p2.ycoord()>y2) y2=p2.ycoord();
  // p3
  if (p3.xcoord()<x1) x1=p3.xcoord();
  if (p3.xcoord()>x2) x2=p3.xcoord(); 
  if (p3.ycoord()<y1) y1=p3.ycoord();
  if (p3.ycoord()>y2) y2=p3.ycoord();
  // p4
  if (p4.xcoord()<x1) x1=p4.xcoord();
  if (p4.xcoord()>x2) x2=p4.xcoord(); 
  if (p4.ycoord()<y1) y1=p4.ycoord();
  if (p4.ycoord()>y2) y2=p4.ycoord();
}

void geowin_Translate(d3_simplex& obj, double dx, double dy)
{
  d3_simplex obnew; 
  int val = GeoWin::get_projection_mode();
  
  switch(val){
   case 0:   { obnew = obj.translate(dx,dy,0); break; }
   case 1:   { obnew = obj.translate(dx,0,dy); break; }
   default: { obnew = obj.translate(0,dx,dy); break; }
  }  
  obj = obnew;
}

void geowin_Rotate(d3_simplex& obj, double dx, double dy,double a){ }

bool geowin_IntersectsBox(const d3_rat_simplex& obj, double x1,double y1,double x2, double y2,bool f)
{ return geowin_IntersectsBox(obj.to_float(),x1,y1,x2,y2,f); }

void geowin_BoundingBox(const d3_rat_simplex& obj, double& x1, double& x2,double& y1, double& y2)
{ geowin_BoundingBox(obj.to_float(),x1,x2,y1,y2); }

void geowin_Translate(d3_rat_simplex& obj, double dx, double dy)
{
 int val = GeoWin::get_projection_mode();
 switch(val){
  case 0: { obj = leda_convert_to(obj.to_float().translate(dx, dy, 0)); break; }
  case 1: { obj = leda_convert_to(obj.to_float().translate(dx, 0, dy)); break; }
  default: { obj = leda_convert_to(obj.to_float().translate(0,dx,dy)); break; }  
 }
}

void geowin_Rotate(d3_rat_simplex& obj, double dx, double dy,double a){ }


void geowin_generate_objects(GeoWin& gw, list<d3_simplex>& L){ }

void geowin_generate_objects(GeoWin& gw, list<d3_rat_simplex>& L){ }


#if (__LEDA__ >= 420)

// d3_triangle / d3_rat_triangle

void d3_triangles_d3(const list<d3_triangle>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_triangle iter;
 forall(iter,L) {
   node v1 = G.new_node(iter.point1());
   node v2 = G.new_node(iter.point2());  
   node v3 = G.new_node(iter.point3());
      
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);
   
   edge e3 = G.new_edge(v1,v3);
   edge e4 = G.new_edge(v3,v1);
   G.set_reversal(e3,e4);
  
   edge e7 = G.new_edge(v2,v3);
   edge e8 = G.new_edge(v3,v2);
   G.set_reversal(e7,e8);
 }
 H.join(G);
}

void d3_rat_triangles_d3(const list<d3_rat_triangle>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_rat_triangle iter;
 forall(iter,L) {
   node v1 = G.new_node(iter.point1().to_float());
   node v2 = G.new_node(iter.point2().to_float());  
   node v3 = G.new_node(iter.point3().to_float());
      
   edge e1 = G.new_edge(v1,v2);
   edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);
   
   edge e3 = G.new_edge(v1,v3);
   edge e4 = G.new_edge(v3,v1);
   G.set_reversal(e3,e4);
   
   edge e7 = G.new_edge(v2,v3);
   edge e8 = G.new_edge(v3,v2);
   G.set_reversal(e7,e8);
 }
 H.join(G);
}

// ---------------------------------------------------------------------------------------
// --------------------------------- generators ------------------------------------------
// ---------------------------------------------------------------------------------------

list<d3_triangle> geowin_conversion(list<d3_rat_point>& L,double wt)
{
 list<d3_triangle> Lf;
 d3_rat_point p;
 d3_point pfl[3];
 int pfl_cnt = 0;
 
 forall(p,L) {
   pfl[pfl_cnt++] = p.to_d3_point().translate(wt,wt,0);
   
   if (pfl_cnt == 3){
    d3_triangle p_tr(pfl[0], pfl[1], pfl[2]);
    Lf.append(p_tr);
    pfl_cnt = 0;
   } 
 }
 return Lf;
}

void geowin_generate_objects(GeoWin& gw, list<d3_triangle>& L)
{
  window& w = gw.get_window();
  int r     = GEOWIN_MARGIN;

  double d = r/w.scale();
  double d1 = d + 16/w.scale();
  double x1 = w.xmin()+d;
  double y1 = w.ymin()+d;
  double x2 = w.xmax()-d;
  double y2 = w.ymax()-d1;

  panel P("D3 Triangles");
  
  double wt = ((x2-x1) > (y2-y1)) ? (y2-y1)/2 : (x2-x1)/2 ;
  int   maxc = (int)wt;
  int   number = 20, number2 =0;
  
  P.int_item("number*100", number2, 0, 19);
  P.int_item("+ n", number, 0, 99);
    
  P.button(" CUBE ", 1);
  P.button(" BALL ", 2);
  P.button(" SQUARE ", 3);

  P.button(" DONE ", 0);

  int but;
  list<d3_rat_point> Lr;

  but = gw.open_panel(P);

  switch (but) {
   case 0: return;
   case 1: random_points_in_cube(3* (number+number2*100),3*maxc/4,Lr); break;
   case 2: random_points_in_ball(3* (number+number2*100),maxc,Lr); break;
   case 3: random_points_in_square(3* (number+number2*100),maxc,Lr); break;
  }

  L=geowin_conversion(Lr,wt);
}

void geowin_generate_objects(GeoWin& gw, list<d3_rat_triangle>& L) 
{ 
  list<d3_triangle> L1;

  geowin_generate_objects(gw, L1);
  d3_triangle t;
  d3_rat_point p1,p2,p3;
  d3_point p;
  rational x,y,z;
  forall(t, L1) {
     p= t.point1();
     x=rational(p.xcoord()); y=rational(p.ycoord());
     z=rational(p.zcoord());
     p1 = d3_rat_point(x,y,z);
     p= t.point2();
     x=rational(p.xcoord()); y=rational(p.ycoord());
     z=rational(p.zcoord());
     p2 = d3_rat_point(x,y,z);
     p= t.point3();
     x=rational(p.xcoord()); y=rational(p.ycoord());
     z=rational(p.zcoord());  
     p3 = d3_rat_point(x,y,z);        
     L.append(d3_rat_triangle(p1,p2,p3)); 
  }
}


ps_file& operator<<(ps_file& F,const d3_triangle& r)
{
   point p1(r.point1().xcoord(), r.point1().ycoord());
   point p2(r.point2().xcoord(), r.point2().ycoord());
   point p3(r.point3().xcoord(), r.point3().ycoord());

   F << segment(p1,p2); F << segment(p1,p3); F << segment(p2,p3);   
   return F;
}

ps_file& operator<<(ps_file& F,const d3_rat_triangle& r)
{
   point p1(r.point1().to_float().xcoord(), r.point1().to_float().ycoord());
   point p2(r.point2().to_float().xcoord(), r.point2().to_float().ycoord());
   point p3(r.point3().to_float().xcoord(), r.point3().to_float().ycoord());

   F << segment(p1,p2); F << segment(p1,p3); F << segment(p2,p3);      
   return F;
}

window& operator << (window& w, const d3_triangle& r)
{
   point p1(r.point1().xcoord(), r.point1().ycoord());
   point p2(r.point2().xcoord(), r.point2().ycoord());
   point p3(r.point3().xcoord(), r.point3().ycoord());

   w << segment(p1,p2); w << segment(p1,p3); w << segment(p2,p3);       
   return w;
}

window& operator << (window& w, const d3_rat_triangle& r)
{
   point p1(r.point1().to_float().xcoord(), r.point1().to_float().ycoord());
   point p2(r.point2().to_float().xcoord(), r.point2().to_float().ycoord());
   point p3(r.point3().to_float().xcoord(), r.point3().to_float().ycoord());

   w << segment(p1,p2); w << segment(p1,p3); w << segment(p2,p3);      
   return w;
}

window& d3_triangle_draw(window& w, const d3_triangle& obj,int val)
{
  val = GeoWin::get_projection_mode();
  point p1,p2,p3;

  switch(val){
   case 0:{
    p1=point(obj.point1().xcoord(), obj.point1().ycoord());
    p2=point(obj.point2().xcoord(), obj.point2().ycoord());
    p3=point(obj.point3().xcoord(), obj.point3().ycoord()); break;
   }   
   case 1:{
    p1=point(obj.point1().xcoord(), obj.point1().zcoord());
    p2=point(obj.point2().xcoord(), obj.point2().zcoord());
    p3=point(obj.point3().xcoord(), obj.point3().zcoord()); break;
   }   
   default:{
    p1=point(obj.point1().ycoord(), obj.point1().zcoord());
    p2=point(obj.point2().ycoord(), obj.point2().zcoord());
    p3=point(obj.point3().ycoord(), obj.point3().zcoord()); break;
   }   
  }
  
  w << segment(p1,p2); w << segment(p1,p3); w << segment(p2,p3); 
  
  return w;
}

window& d3_rat_triangle_draw(window& w, const d3_rat_triangle& obj,int val)
{
  return d3_triangle_draw(w,obj.to_float(),val);
}

window& operator >> (window& w, d3_triangle& obj)
{
   d3_point p1,p2,p3;
   point p;
   w >> p;
   double x,y;
   
   p1 = d3_point(p.xcoord(),p.ycoord(),0);
   w.read_mouse_seg(p.xcoord(),p.ycoord(),x,y);
   w << segment(p.xcoord(),p.ycoord(),x,y);
   p2 = d3_point(x,y,0);
   w.read_mouse_seg(p2.xcoord(),p2.ycoord(),x,y);
   p3 = d3_point(x,y,0);
   
   obj = d3_triangle(p1,p2,p3);
   return w;
}

window& operator >> (window& w, d3_rat_triangle& obj)
{
   d3_triangle t;
   w >> t;
   obj = leda_convert_to(t);
   return w;
}

bool geowin_IntersectsBox(const d3_triangle& obj, double x1,double y1,double x2, double y2,bool f)
{
  int val = GeoWin::get_projection_mode();
  point p1,p2,p3;

  switch(val){
   case 0:{
    p1=point(obj.point1().xcoord(), obj.point1().ycoord());
    p2=point(obj.point2().xcoord(), obj.point2().ycoord());
    p3=point(obj.point3().xcoord(), obj.point3().ycoord()); break;
   }   
   case 1:{
    p1=point(obj.point1().xcoord(), obj.point1().zcoord());
    p2=point(obj.point2().xcoord(), obj.point2().zcoord());
    p3=point(obj.point3().xcoord(), obj.point3().zcoord()); break;
   }   
   default:{
    p1=point(obj.point1().ycoord(), obj.point1().zcoord());
    p2=point(obj.point2().ycoord(), obj.point2().zcoord());
    p3=point(obj.point3().ycoord(), obj.point3().zcoord()); break;
   }   
  }
    
  rectangle R(x1,y1,x2,y2);
  list<point> IL;
  IL = R.intersection(segment(p1,p2));
  if (! IL.empty()) return true;
  IL = R.intersection(segment(p1,p3));
  if (! IL.empty()) return true;  
  IL = R.intersection(segment(p2,p3));
  if (! IL.empty()) return true;
  return false;
}

void geowin_BoundingBox(const d3_triangle& obj, double& x1, double& x2,double& y1, double& y2)
{  
  int val = GeoWin::get_projection_mode();
  point p1,p2,p3,p4;
  switch(val){
   case 0:{
    p1=point(obj.point1().xcoord(), obj.point1().ycoord());
    p2=point(obj.point2().xcoord(), obj.point2().ycoord());
    p3=point(obj.point3().xcoord(), obj.point3().ycoord()); break;
   }   
   case 1:{
    p1=point(obj.point1().xcoord(), obj.point1().zcoord());
    p2=point(obj.point2().xcoord(), obj.point2().zcoord());
    p3=point(obj.point3().xcoord(), obj.point3().zcoord()); break;
   }   
   default:{
    p1=point(obj.point1().ycoord(), obj.point1().zcoord());
    p2=point(obj.point2().ycoord(), obj.point2().zcoord());
    p3=point(obj.point3().ycoord(), obj.point3().zcoord()); break;
   }   
  } 
  
  x1=p1.xcoord(); x2=p1.xcoord(); y1=p1.ycoord(); y2=p1.ycoord();
  // p2
  if (p2.xcoord()<x1) x1=p2.xcoord();
  if (p2.xcoord()>x2) x2=p2.xcoord(); 
  if (p2.ycoord()<y1) y1=p2.ycoord();
  if (p2.ycoord()>y2) y2=p2.ycoord();
  // p3
  if (p3.xcoord()<x1) x1=p3.xcoord();
  if (p3.xcoord()>x2) x2=p3.xcoord(); 
  if (p3.ycoord()<y1) y1=p3.ycoord();
  if (p3.ycoord()>y2) y2=p3.ycoord();
}

void geowin_Translate(d3_triangle& obj, double dx, double dy)
{
  d3_triangle obnew; 
  int val = GeoWin::get_projection_mode();
  
  switch(val){
   case 0:   { obnew = obj.translate(dx,dy,0); break; }
   case 1:   { obnew = obj.translate(dx,0,dy); break; }
   default: { obnew = obj.translate(0,dx,dy); break; }
  }  
  obj = obnew;
}

void geowin_Rotate(d3_triangle& obj, double dx, double dy,double a){ }

bool geowin_IntersectsBox(const d3_rat_triangle& obj, double x1,double y1,double x2, double y2,bool f)
{ return geowin_IntersectsBox(obj.to_float(),x1,y1,x2,y2,f); }

void geowin_BoundingBox(const d3_rat_triangle& obj, double& x1, double& x2,double& y1, double& y2)
{ geowin_BoundingBox(obj.to_float(),x1,x2,y1,y2); }

void geowin_Translate(d3_rat_triangle& obj, double dx, double dy)
{
 int val = GeoWin::get_projection_mode();
 switch(val){
  case 0: { obj = leda_convert_to(obj.to_float().translate(dx, dy, 0)); break; }
  case 1: { obj = leda_convert_to(obj.to_float().translate(dx, 0, dy)); break; }
  default: { obj = leda_convert_to(obj.to_float().translate(0,dx,dy)); break; }  
 }
}

void geowin_Rotate(d3_rat_triangle& obj, double dx, double dy,double a){ }




// ------------------------------------------------------------------------------
// New  - d3 polygons
// ------------------------------------------------------------------------------

/*
#if !defined(__SUNPRO_CC) 
void d3_polygons_d3(const list<d3_polygon>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_polygon iter;
 list<d3_segment> LS;
 forall(iter,L) {
   LS = iter.segments();
   d3_segment siter;
   forall(siter,LS) {
    node v1 = G.new_node(siter.source());
    node v2 = G.new_node(siter.target());   
    edge e1 = G.new_edge(v1,v2);
    edge e2 = G.new_edge(v2,v1);
    G.set_reversal(e1,e2);
   } 
 }
 H.join(G);
}

void d3_rat_polygons_d3(const list<d3_rat_polygon>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 d3_rat_polygon iter;
 list<d3_rat_segment> LS;
 forall(iter,L) {
   LS = iter.segments();
   d3_rat_segment siter;
   forall(siter,LS) {
    node v1 = G.new_node(siter.source().to_float());
    node v2 = G.new_node(siter.target().to_float());   
    edge e1 = G.new_edge(v1,v2);
    edge e2 = G.new_edge(v2,v1);
    G.set_reversal(e1,e2);
   } 
 }
 H.join(G);
}

ps_file& operator<<(ps_file& F,const d3_polygon& p)
{
   return F;
}

ps_file& operator<<(ps_file& F,const d3_rat_polygon& p)
{     
   return F;
}

window& operator << (window& w, const d3_polygon& p)
{
   list<d3_segment> LS = p.segments();
   d3_segment siter;
   forall(siter, LS){
     w << siter.project_xy();
   }
   return w;
}

window& operator << (window& w, const d3_rat_polygon& p)
{
   list<d3_rat_segment> LS = p.segments();
   d3_rat_segment siter;
   forall(siter, LS){
     w << siter.project_xy().to_float();
   } 
   return w;
}

window& operator >> (window& w, d3_polygon& obj)
{
   polygon p;
   w >> p;
   list<d3_point> LP;
   
   list<point> Lh = p.vertices();
   point pi;
   forall(pi,Lh) LP.append(d3_point(pi.xcoord(), pi.ycoord(), 0));
     
   obj = d3_polygon(LP);
   return w;
}

window& operator >> (window& w, d3_rat_polygon& obj)
{
   d3_polygon p;
   w >> p;
   obj = leda_convert_to(p);
   return w;
}


bool geowin_IntersectsBox(const d3_polygon& obj, double x1,double y1,double x2, double y2,bool f)
{
  list<d3_point> Lp = obj.vertices();
  list<point> Lh;
  d3_point pi;
  forall(pi,Lp) Lh.append(pi.project_xy());
  
  polygon p(Lh);

  return geowin_IntersectsBox(p, x1,y1,x2,y2,f);
}

void geowin_BoundingBox(const d3_polygon& obj, double& x1, double& x2,double& y1, double& y2)
{  
  list<d3_point> Lp = obj.vertices();
  list<point> Lh;
  d3_point pi;
  forall(pi,Lp) Lh.append(pi.project_xy());
  
  polygon p(Lh);

  geowin_BoundingBox(p, x1,x2,y1,y2);
}

void geowin_Translate(d3_polygon& obj, double dx, double dy)
{
  d3_polygon obnew; 
  int val = GeoWin::get_projection_mode();
  
  switch(val){
   case 0:   { obnew = obj.translate(dx,dy,0); break; }
   case 1:   { obnew = obj.translate(dx,0,dy); break; }
   default: { obnew = obj.translate(0,dx,dy); break; }
  }  
  obj = obnew;
}

void geowin_Rotate(d3_polygon& obj, double dx, double dy,double a){ }

bool geowin_IntersectsBox(const d3_rat_polygon& obj, double x1,double y1,double x2, double y2,bool f)
{ return geowin_IntersectsBox(obj.to_float(),x1,y1,x2,y2,f); }

void geowin_BoundingBox(const d3_rat_polygon& obj, double& x1, double& x2,double& y1, double& y2)
{ geowin_BoundingBox(obj.to_float(),x1,x2,y1,y2); }

void geowin_Translate(d3_rat_polygon& obj, double dx, double dy)
{
 int val = GeoWin::get_projection_mode();
 switch(val){
  case 0: { obj = leda_convert_to(obj.to_float().translate(dx, dy, 0)); break; }
  case 1: { obj = leda_convert_to(obj.to_float().translate(dx, 0, dy)); break; }
  default: { obj = leda_convert_to(obj.to_float().translate(0,dx,dy)); break; }  
 }
}

void geowin_Rotate(d3_rat_polygon& obj, double dx, double dy,double a){ }


void geowin_generate_objects(GeoWin& gw, list<d3_polygon>& L){ }

void geowin_generate_objects(GeoWin& gw, list<d3_rat_polygon>& L){ }

#endif
*/


#endif


// d3_sphere / d3_rat_sphere
// d3 functions
d3_sphere convert_to_float(const d3_rat_sphere& sp)
{
  return d3_sphere(sp.point1().to_float(), sp.point2().to_float(), sp.point3().to_float(), sp.point4().to_float());
}

void d3_spheres_d3(const list<d3_sphere>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
 GRAPH<d3_point,int> G;
 GRAPH<d3_point,int> G2;
 d3_sphere iter;
 list<d3_point> loc;
 list<d3_rat_point> loc_rat;
 //random_points_on_sphere(200,100,loc_rat);
 generate_points_on_sphere(d3_sphere(d3_point(100,100,100),d3_point(-100,100,100),d3_point(100,-100,100),d3_point(-100,-100,-100)), loc);
 
 //d3_rat_point pd3;
 //forall(pd3,loc_rat) loc.append(pd3.to_float());
 //cout << "chull !\n";
 CONVEX_HULL(loc,G);
 
 vector v;
 d3_point orig(0,0,0);
 double scale;
 node piter;
 
 forall(iter,L) {
    d3_point ct = iter.center();
    double rd = iter.radius();
    //move vector origin->center
    v = ct - orig;
    scale = rd/173.2;
    G2.clear();
    G2=G;
    
    forall_nodes(piter,G2) {
      d3_point hp(G2[piter].xcoord()*scale, G2[piter].ycoord()*scale, G2[piter].zcoord()*scale);
      G2[piter]= hp + v;
    }

    H.join(G2);
    G2.clear();
 }
}

void d3_rat_spheres_d3(const list<d3_rat_sphere>& L,d3_window& W, GRAPH<d3_point,int>& H)
{
  list<d3_sphere> Lf;
  d3_rat_sphere iter;
  forall(iter,L) Lf.append(convert_to_float(iter));
  d3_spheres_d3(Lf,W,H);
}

ps_file& operator<<(ps_file& F,const d3_sphere& S)
{
   d3_point pc = S.center();
   point p1(pc.xcoord(),pc.ycoord()); //center
   double rd = S.radius();
   point p2(pc.xcoord()-rd,pc.ycoord());
   
   F << circle(p1,p2);  
   return F;
}

ps_file& operator<<(ps_file& F,const d3_rat_sphere& S)
{
   F << convert_to_float(S);
   return F;
}


window& operator << (window& w, const d3_sphere& S)
{
   d3_point pc = S.center();
   point p1(pc.xcoord(),pc.ycoord()); //center
   double rd = S.radius();
   point p2(pc.xcoord()-rd,pc.ycoord());
   
   w << circle(p1,p2);  
   return w;
}

window& operator << (window& w, const d3_rat_sphere& S)
{
   w << convert_to_float(S);
   return w;
}

window& d3_sphere_draw(window& w, const d3_sphere& obj,int val)
{
  val = GeoWin::get_projection_mode();
  d3_point pc=obj.center();
  double rd=obj.radius();
  point p1,p2;

  switch(val){
   case 0:{
    p1=point(pc.xcoord(),pc.ycoord()); p2=point(pc.xcoord()-rd,pc.ycoord()); break;
   }   
   case 1:{
    p1=point(pc.xcoord(),pc.zcoord()); p2=point(pc.xcoord()-rd,pc.zcoord()); break;
   }   
   default:{
    p1=point(pc.ycoord(),pc.zcoord()); p2=point(pc.ycoord()-rd,pc.zcoord()); break;
   }   
  }
  
  w << circle(p1,p2);
  
  return w;
}

window& d3_rat_sphere_draw(window& w, const d3_rat_sphere& obj,int val)
{
  return d3_sphere_draw(w,obj.to_float(),val);
}


window& operator >> (window& w, d3_sphere& obj)
{
   circle ch;
   w >> ch;
   double r  = ch.radius();
   point   cp = ch.center();
   double xc = cp.xcoord(), yc = cp.ycoord();
   d3_point p1(xc+r,yc,0), p2(xc,yc+r,0);
   d3_point p3(xc,yc,r), p4(xc,yc-r,0);
   
   obj= d3_sphere(p1,p2,p3,p4);
   return w;
}

window& operator >> (window& w, d3_rat_sphere& obj)
{
   d3_sphere h;
   w >> h;
   obj = leda_convert_to(h);
   return w;
}

bool geowin_IntersectsBox(const d3_sphere& obj, double x1,double y1,double x2, double y2,bool f)
{
   int val = GeoWin::get_projection_mode();
   d3_point pc = obj.center();
   point p1,p2;
   double rd = obj.radius();
 
   switch(val){
   case 0:{
    p1=point(pc.xcoord(),pc.ycoord()); p2=point(pc.xcoord()-rd,pc.ycoord()); break;
   }   
   case 1:{
    p1=point(pc.xcoord(),pc.zcoord()); p2=point(pc.xcoord()-rd,pc.zcoord()); break;
   }   
   default:{
    p1=point(pc.ycoord(),pc.zcoord()); p2=point(pc.ycoord()-rd,pc.zcoord()); break;
   }   
   }

   circle Cf(p1,p2);  
   return geowin_IntersectsBox(Cf,x1,y1,x2,y2,f);
}

void geowin_BoundingBox(const d3_sphere& obj, double& x1, double& x2,double& y1, double& y2)
{  
   int val = GeoWin::get_projection_mode();
   d3_point pc = obj.center();
   point p1,p2;
   double rd = obj.radius();

   switch(val){
   case 0:{
    p1=point(pc.xcoord(),pc.ycoord()); p2=point(pc.xcoord()-rd,pc.ycoord()); break;
   }   
   case 1:{
    p1=point(pc.xcoord(),pc.zcoord()); p2=point(pc.xcoord()-rd,pc.zcoord()); break;
   }   
   default:{
    p1=point(pc.ycoord(),pc.zcoord()); p2=point(pc.ycoord()-rd,pc.zcoord()); break;
   }   
   }
   
   circle Cf(p1,p2);  
   geowin_BoundingBox(Cf,x1,y1,x2,y2);
}

void geowin_Translate(d3_sphere& obj, double dx, double dy)
{
  d3_sphere obnew; 
  int val = GeoWin::get_projection_mode();
  
  switch(val){
   case 0:   { obnew = obj.translate(dx,dy,0); break; }
   case 1:   { obnew = obj.translate(dx,0,dy); break; }
   default: { obnew = obj.translate(0,dx,dy); break; }
  }  
  obj = obnew;
}

void geowin_Rotate(d3_sphere& obj, double dx, double dy,double a){ }

bool geowin_IntersectsBox(const d3_rat_sphere& obj, double x1,double y1,double x2, double y2,bool f)
{ return geowin_IntersectsBox(convert_to_float(obj),x1,y1,x2,y2,f); }

void geowin_BoundingBox(const d3_rat_sphere& obj, double& x1, double& x2,double& y1, double& y2)
{ geowin_BoundingBox(convert_to_float(obj),x1,x2,y1,y2); }

void geowin_Translate(d3_rat_sphere& obj, double dx, double dy)
{
 int val = GeoWin::get_projection_mode();
 switch(val){
  case 0: { obj = leda_convert_to(obj.to_float().translate(dx, dy, 0)); break; }
  case 1: { obj = leda_convert_to(obj.to_float().translate(dx, 0, dy)); break; }
  default: { obj = leda_convert_to(obj.to_float().translate(0,dx,dy)); break; }  
 }
}

void geowin_Rotate(d3_rat_sphere& obj, double dx, double dy,double a){ }


void geowin_generate_objects(GeoWin& gw, list<d3_sphere>& L){ }

void geowin_generate_objects(GeoWin& gw, list<d3_rat_sphere>& L){ }

#if LEDA_ROOT_INCL_ID == 400955
#undef LEDA_ROOT_INCL_ID
#include <LEDA/UNDEFINE_NAMES.h>
#endif

#endif

