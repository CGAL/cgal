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
// file          : src/GeoWin/geo_ledaobj.c
// package       : GeoWin (1.2.2)
// revision      : 1.2.2
// revision_date : 30 January 2000 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ============================================================================

#include <LEDA/geowin_init.h>

GEOWIN_BEGIN_NAMESPACE

void geowin_Translate(point& obj, double dx, double dy)
{ 
  point obnew = obj.translate(dx,dy);
  obj = obnew;
}

void geowin_Translate(segment& obj, double dx, double dy)
{ 
  segment obnew = obj.translate(dx,dy);
  obj = obnew;
}

void geowin_Translatepoint_seg(segment& obj, double dx, double dy, int pnr)
{ 
  segment obnew;
  if (pnr==0) obnew = segment(obj.source().translate(dx,dy),obj.target()); 
  else obnew = segment(obj.source(),obj.target().translate(dx,dy)); 
  obj = obnew;
}

void geowin_Translate(ray& obj, double dx, double dy)
{ 
  ray obnew = obj.translate(dx,dy);
  obj = obnew;
}

void geowin_Translatepoint_ray(ray& obj, double dx, double dy, int pnr)
{ 
  ray obnew;
  if (pnr==0) obnew = ray(obj.point1().translate(dx,dy),obj.point2()); 
  else obnew = ray(obj.point1(),obj.point2().translate(dx,dy)); 
  obj = obnew;
}

void geowin_Translate(line& obj, double dx, double dy)
{ 
  line obnew = obj.translate(dx,dy);
  obj = obnew;
}

void geowin_Translatepoint_line(line& obj, double dx, double dy, int pnr)
{ 
  line obnew;
  if (pnr==0) obnew = line(obj.point1().translate(dx,dy),obj.point2()); 
  else obnew = line(obj.point1(),obj.point2().translate(dx,dy)); 
  obj = obnew;
}

void geowin_Translate(circle& obj, double dx, double dy)
{ 
  circle obnew = obj.translate(dx,dy);
  obj = obnew;
}

void geowin_Translatepoint_circle(circle& obj, double dx, double dy, int pnr)
{ 
  circle obnew;
  if (pnr==0) obnew = circle(obj.point1().translate(dx,dy),obj.point2(),obj.point3()); 
  else {
   if (pnr==1) obnew = circle(obj.point1(),obj.point2().translate(dx,dy),obj.point3()); 
   else
    obnew = circle(obj.point1(),obj.point2(),obj.point3().translate(dx,dy)); 
  }
  obj = obnew;
}


void geowin_Translate(polygon& obj, double dx, double dy)
{ 
  polygon obnew = obj.translate(dx,dy);
  obj = obnew;
}

void geowin_Translatepoint_poly(polygon& obj, double dx, double dy, int pnr)
{ 
  list<point> PL = obj.vertices();
  if (pnr >= PL.size()) return;
  // translate a polygon vertex 
  list_item it = PL.get_item(pnr);
  //cout << PL[it] << " " << pnr << "\n"; cout.flush();
  PL[it]=PL[it].translate(dx,dy);
  obj = polygon(PL);
}

void geowin_Translate(gen_polygon& obj, double dx, double dy)
{ 
  gen_polygon obnew = obj.translate(dx,dy);
  obj = obnew;
}

void geowin_Translatepoint_gpoly(gen_polygon& obj, double dx, double dy, int pnr)
{ 
  list<point> PL = obj.vertices();
  if (pnr >= PL.size()) return;
  // translate a polygon vertex 
  list_item it = PL.get_item(pnr);
  //cout << PL[it] << " " << pnr << "\n"; cout.flush();
  PL[it]=PL[it].translate(dx,dy);
  obj = gen_polygon(PL);  
}

void geowin_Translate(d3_point& obj, double dx, double dy)
{ 
  int dm = GeoWin::get_projection_mode();
  d3_point obnew;
  
  switch(dm){
   case 0: { obnew = obj.translate(dx,dy,0); break; } //xy
   case 1: { obnew = obj.translate(dx,0,dy); break; } //xz
   default: { obnew = obj.translate(0,dx,dy); break; } //yz
  }
  
  obj = obnew;
}

void geowin_Translate(rectangle& obj, double dx, double dy)
{ 
  rectangle obnew = obj.translate(dx,dy);
  obj = obnew;
}

#if (__LEDA__ >= 420)
void geowin_Translate(triangle& obj, double dx, double dy)
{ 
  triangle obnew = obj.translate(dx,dy);
  obj = obnew;
}
#endif


void geowin_Rotate(point& obj, double x, double y, double a)
{ 
 point objnew = obj.rotate(point(x,y), a); 
 obj = objnew;
}

void geowin_Rotate(segment& obj, double x, double y, double a)
{ 
 segment objnew = obj.rotate(point(x,y), a); 
 obj = objnew;
}

void geowin_Rotate(ray& obj, double x, double y, double a)
{ 
 ray objnew = obj.rotate(point(x,y), a); 
 obj = objnew;
}

void geowin_Rotate(line& obj, double x, double y, double a)
{ 
 line objnew = obj.rotate(point(x,y), a); 
 obj = objnew;
}

void geowin_Rotate(circle& obj, double x, double y, double a)
{ 
 circle objnew = obj.rotate(point(x,y), a); 
 obj = objnew;
}

void geowin_Rotate(polygon& obj, double x, double y, double a)
{ 
 polygon objnew = obj.rotate(point(x,y), a); 
 obj = objnew;
}

void geowin_Rotate(gen_polygon& obj, double x, double y, double a)
{ 
 gen_polygon objnew = obj.rotate(point(x,y), a); 
 obj = objnew;
}

void geowin_Rotate(d3_point& obj, double x, double y, double a)
{ return ; }

void geowin_Rotate(rectangle& obj, double x, double y, double a)
{ return; } 

#if (__LEDA__ >= 420)
void geowin_Rotate(triangle& obj, double x, double y, double a)
{ 
 triangle objnew = obj.rotate(point(x,y), a); 
 obj = objnew;
}
#endif


// IntersectsBox and BoundingBox for leda-types

// ************************************************************************
// ******************** begin intersects **********************************

bool geowin_IntersectsBox(const point& p, double x1, double y1,
		       double x2, double y2,bool filled)
{ 
  double x = p.xcoord(), y = p.ycoord();
  
  if ( x < x1 || x > x2 ) return false;
  if ( y < y1 || y > y2 ) return false;
  
  return true; 
}

bool geowin_IntersectsBox(const segment& s, double x1, double y1,
		       double x2, double y2,bool filled)
{
  segment s1(x1, y1, x2, y1);
  segment s2(x1, y2, x2, y2);
  segment s3(x1, y1, x1, y2);
  segment s4(x2, y1, x2, y2);
  point p;
  
  return 
    ( geowin_IntersectsBox(s.start(), x1, y1, x2, y2,filled) ||
      geowin_IntersectsBox(s.end(), x1, y1, x2, y2,filled)   ||
      s1.intersection(s,p) || s2.intersection(s,p) ||
      s3.intersection(s,p) || s4.intersection(s,p) );
}


bool geowin_IntersectsBox(const ray& r, double x1, double y1,
		       double x2, double y2,bool filled )
{
  segment s1(x1, y1, x2, y1);
  segment s2(x1, y2, x2, y2);
  segment s3(x1, y1, x1, y2);
  segment s4(x2, y1, x2, y2);
  point p;
  
  return 
    (  r.intersection(s1,p) || r.intersection(s2,p) ||
       r.intersection(s3,p) || r.intersection(s4,p)  );
}


bool geowin_IntersectsBox(const line& l, double x1, double y1,
		       double x2, double y2 ,bool filled)
{ 
  segment s1(x1, y1, x2, y1);
  segment s2(x1, y2, x2, y2);
  segment s3(x1, y1, x1, y2);
  segment s4(x2, y1, x2, y2);
  point p;
  
  return 
    (  l.intersection(s1,p) || l.intersection(s2,p) ||
       l.intersection(s3,p) || l.intersection(s4,p)  );
}

bool geowin_IntersectsBox(const circle& c, double x1, double y1,
		       double x2, double y2 ,bool filled)
{
  if (!c.orientation())
    return geowin_IntersectsBox(segment(c.point1(), c.point2()), x1, y1, x2, y2, filled );
  
  segment s1(x1, y1, x2, y1);
  segment s2(x1, y2, x2, y2);
  segment s3(x1, y1, x1, y2);
  segment s4(x2, y1, x2, y2);
  
  if ( !c.intersection(s1).empty() ) return true;
  if ( !c.intersection(s2).empty() ) return true;
  if ( !c.intersection(s3).empty() ) return true;
  if ( !c.intersection(s4).empty() ) return true;
  
  if( c.inside(point(x1, y1)) ) return filled;
  return geowin_IntersectsBox(c.center(), x1, y1, x2, y2, filled);
}

bool geowin_IntersectsBox(const polygon& p,  double x1, double y1,
		       double x2, double y2 , bool filled)
{
  if (p.size() == 0) return false;
  if (p.size() == 1) geowin_IntersectsBox(p.vertices().head(), x1, y1, x2, y2 ,filled);

  segment s1(x1, y1, x2, y1);
  segment s2(x1, y2, x2, y2);
  segment s3(x1, y1, x1, y2);
  segment s4(x2, y1, x2, y2);
  
  int ori = p.orientation();

  if ( !p.intersection(s1).empty() ) return true;
  if ( !p.intersection(s2).empty() ) return true;
  if ( !p.intersection(s3).empty() ) return true;
  if ( !p.intersection(s4).empty() ) return true;
  
  if (ori == -1) {
   if( p.outside(point(x1, y1)) ) return filled;  
  }
  else {
   if( p.inside(point(x1, y1)) ) return filled;
  }
  
  return geowin_IntersectsBox(p.vertices().head(), x1, y1, x2, y2, filled);
}

bool geowin_IntersectsBox(const gen_polygon& p,  double x1, double y1,
		       double x2, double y2 , bool filled)
{
  segment s1(x1, y1, x2, y1);
  segment s2(x1, y2, x2, y2);
  segment s3(x1, y1, x1, y2);
  segment s4(x2, y1, x2, y2);
  
  if ( !p.intersection(s1).empty() ) return true;
  if ( !p.intersection(s2).empty() ) return true;
  if ( !p.intersection(s3).empty() ) return true;
  if ( !p.intersection(s4).empty() ) return true;
  
  if( p.inside(point(x1, y1)) ) return filled;

  list<polygon>  pli= p.polygons();
  polygon pl;

  forall(pl,pli) {
   if (geowin_IntersectsBox(pl.vertices().head(), x1, y1, x2, y2, filled)) return true;
  }
  return false;
}

#if (__LEDA__ >= 420)
bool geowin_IntersectsBox(const triangle& t,  double x1, double y1,
		       double x2, double y2 , bool filled)
{
  list<point> LP;
  if (left_turn(t.point1(),t.point2(),t.point3()))
   { LP.append(t.point1()); LP.append(t.point2()); LP.append(t.point3()); }
  else
   { LP.push(t.point1()); LP.push(t.point2()); LP.push(t.point3()); } 
  polygon pol(LP);
  return geowin_IntersectsBox(pol,x1,y1,x2,y2,filled);
}
#endif

bool geowin_IntersectsBox(const rectangle& obj, double x1,double y1,double x2, double y2,bool f)
{
 point p1(x1,y1),p2(x2,y2);

 list<rectangle> lr=obj.intersection(rectangle(p1,p2));
 if (lr.size()>0 ) return true; else return false;
}

bool geowin_IntersectsBox(const d3_point& p,double x1,double y1,double x2, double y2,bool f)
{ 
  int dm = GeoWin::get_projection_mode();
  double xw,yw;
  
  switch(dm){
    case 0:  { xw=p.xcoord(); yw=p.ycoord(); break; } //xy
    case 1:  { xw=p.xcoord(); yw=p.zcoord(); break; } //xz
    default:{ xw=p.ycoord(); yw=p.zcoord(); break; } //yz
  }
  
  if (x1<=xw && x2>=xw && y1<=yw && y2>=yw) return true; else return false;
}


// ********************* end  intersects **********************************
// ************************************************************************


// ************************************************************************
// ******************** begin BB *****************************************

void geowin_BoundingBox(const point& p, double& x1, double& x2, double& y1, 
		double& y2)
{
   x1 = p.xcoord(); x2 = p.xcoord();
   y1 = p.ycoord(); y2 = p.ycoord();
}

void geowin_BoundingBox(const segment& s, double& x1, double& x2, 
		double& y1, double& y2)
{
  if (s.xcoord1()<s.xcoord2()) { x1=s.xcoord1(); x2=s.xcoord2(); }
  else { x2=s.xcoord1(); x1=s.xcoord2(); }

  if (s.ycoord1()<s.ycoord2()) { y1=s.ycoord1(); y2=s.ycoord2(); }
  else { y2=s.ycoord1(); y1=s.ycoord2(); }

}

void geowin_BoundingBox(const ray& r, double& x1, double& x2,
		double& y1, double& y2)
{ 
  geowin_BoundingBox(segment(r.source(),r.point2()), x1, x2, y1, y2);
}

void geowin_BoundingBox(const line& l, double& x1, double& x2, double& y1, double& y2)
{ 
  geowin_BoundingBox(segment(l.point1(),l.point2()), x1, x2, y1, y2);
}

void geowin_BoundingBox(const circle& c, double& x1, double& x2, 
		double& y1, double& y2)
{
  point mp = c.center();
  double cx=mp.xcoord();
  double cy=mp.ycoord();
  double r = c.radius();
  
  x1= cx-r; x2=cx+r; y1=cy-r; y2=cy+r;
}
 

void geowin_BoundingBox(const polygon& P, double& x1, double& x2,
		double& y1, double& y2)
{
  point p;
  const list<point>&  L = P.vertices(); 

  if (L.empty()) return;

  x1= (L.head()).xcoord(); x2=x1;
  y1= (L.head()).ycoord(); y2=y1;
  double xakt,yakt;

  forall(p, L){
   xakt=p.xcoord(); yakt=p.ycoord();
   if (xakt<x1) x1=xakt;
   if (xakt>x2) x2=xakt;
   if (yakt<y1) y1=yakt;
   if (yakt>y2) y2=yakt;
  }
}

void geowin_BoundingBox(const gen_polygon& P, double& x1, double& x2,
		double& y1, double& y2)
{
  point p;
  const list<point>&  L = P.vertices(); // polygons with 0 points ??
  x1= (L.head()).xcoord(); x2=x1;
  y1= (L.head()).ycoord(); y2=y1;
  double xakt,yakt;

  forall(p, L){
   xakt=p.xcoord(); yakt=p.ycoord();
   if (xakt<x1) x1=xakt;
   if (xakt>x2) x2=xakt;
   if (yakt<y1) y1=yakt;
   if (yakt>y2) y2=yakt;
  }
}

#if (__LEDA__ >= 420)
void geowin_BoundingBox(const triangle& t, double& x1, double& x2,
		double& y1, double& y2)
{
  point p;
  list<point> L;
  L.append(t.point1()); L.append(t.point2()); L.append(t.point3());
  x1= (L.head()).xcoord(); x2=x1;
  y1= (L.head()).ycoord(); y2=y1;
  double xakt,yakt;

  forall(p, L){
   xakt=p.xcoord(); yakt=p.ycoord();
   if (xakt<x1) x1=xakt;
   if (xakt>x2) x2=xakt;
   if (yakt<y1) y1=yakt;
   if (yakt>y2) y2=yakt;
  }
}
#endif

void geowin_BoundingBox(const rectangle& obj, double& x1, double& x2,double& y1, double& y2)
{
  x1=obj.xmin(); x2=obj.xmax(); y1=obj.ymin(); y2=obj.ymax(); 
}

void geowin_BoundingBox(const d3_point& p, double& x1, double& x2, double& y1, double& y2)
{ 
  int dm = GeoWin::get_projection_mode();
  switch(dm){
    case 0:   { x1=p.xcoord(); x2=p.xcoord(); y1=p.ycoord(); y2=p.ycoord(); return; }
    case 1:   { x1=p.xcoord(); x2=p.xcoord(); y1=p.zcoord(); y2=p.zcoord(); return; }
    default: { x1=p.xcoord(); x2=p.xcoord(); y1=p.zcoord(); y2=p.zcoord(); return; }
  }
}


// ********************* end BB   *****************************************
// ************************************************************************





#if !defined(NO_RAT_ALGORITHMS)

void geowin_Translate(rat_point& p, double dx, double dy)
{ p = rat_point(p.to_point().translate(dx, dy), 20); }

void geowin_Translate(rat_segment& s, double dx, double dy)
{ s = rat_segment(s.to_segment().translate(dx, dy), 20); }

void geowin_Translatepoint_rat_seg(rat_segment& s, double dx, double dy,int pnr)
{ segment sf = s.to_segment(); 
  geowin_Translatepoint_seg(sf,dx,dy,pnr);
  s = rat_segment(sf, 20); 
}

void geowin_Translate(rat_ray& r, double dx, double dy)
{ r = rat_ray(r.to_ray().translate(dx, dy), 20); }

void geowin_Translatepoint_rat_ray(rat_ray& r, double dx, double dy,int pnr)
{ ray rf = r.to_ray(); 
  geowin_Translatepoint_ray(rf,dx,dy,pnr);
  r = rat_ray(rf, 20); 
}

void geowin_Translate(rat_line& l, double dx, double dy)
{ l = rat_line(l.to_line().translate(dx, dy), 20); }

void geowin_Translatepoint_rat_line(rat_line& l, double dx, double dy,int pnr)
{ line lf = l.to_line(); 
  geowin_Translatepoint_line(lf,dx,dy,pnr);
  l = rat_line(lf, 20); 
}

void geowin_Translate(rat_circle& ci, double dx, double dy)
{ ci = rat_circle(ci.to_circle().translate(dx, dy), 20); }

void geowin_Translatepoint_rat_circle(rat_circle& c, double dx, double dy,int pnr)
{ circle cf = c.to_circle(); 
  geowin_Translatepoint_circle(cf,dx,dy,pnr);
  c = rat_circle(cf, 20); 
}

void geowin_Translate(rat_polygon& p, double dx, double dy)
{ p = rat_polygon(p.to_polygon().translate(dx, dy), 20); }

void geowin_Translatepoint_rat_poly(rat_polygon& p, double dx, double dy, int pnr)
{ polygon pf = p.to_polygon(); 
  geowin_Translatepoint_poly(pf,dx,dy,pnr);
  p = rat_polygon(pf, 20); 
}

void geowin_Translate(rat_gen_polygon& p, double dx, double dy)
{ p = rat_gen_polygon(p.to_gen_polygon().translate(dx, dy), 20); }

void geowin_Translatepoint_rat_gpoly(rat_gen_polygon& p, double dx, double dy,int pnr)
{ gen_polygon pf = p.to_gen_polygon(); 
  geowin_Translatepoint_gpoly(pf,dx,dy,pnr);
  p = rat_gen_polygon(pf, 20); 
}

void geowin_Translate(d3_rat_point& p, double dx, double dy)
{ 
  int dm = GeoWin::get_projection_mode();
  d3_point phelp;
  
  switch(dm){
    case 0:  { phelp= p.to_d3_point().translate(dx,dy,0); break; }
    case 1:  { phelp= p.to_d3_point().translate(dx,0,dy); break; }
    default:{ phelp= p.to_d3_point().translate(0,dx,dy); break; }
  }
  
  double x=phelp.xcoord(), y=phelp.ycoord(), z=phelp.zcoord();
  p = d3_rat_point(integer(x*100000),integer(y*100000), integer(z*100000),                                                                  100000);
}

#if (__LEDA__ >= 420)
void geowin_Translate(rat_triangle& obj, double dx, double dy)
{
  triangle t = obj.to_float().translate(dx,dy);
  obj = rat_triangle(rat_point(t.point1(),20), rat_point(t.point2(),20), rat_point(t.point3(),20));
}
#endif

void geowin_Translate(rat_rectangle& obj, double dx, double dy)
{  
  obj = rat_rectangle(obj.to_float().translate(dx, dy), 20);
}


void geowin_Rotate(rat_point& p, double x, double y,double a)
{ p = rat_point(p.to_point().rotate(point(x, y), a), 20); }

void geowin_Rotate(rat_segment& s, double x, double y,double a)
{ s = rat_segment(s.to_segment().rotate(point(x, y), a), 20); }

void geowin_Rotate(rat_ray& r, double x, double y,double a)
{ r = rat_ray(r.to_ray().rotate(point(x, y), a), 20); }

void geowin_Rotate(rat_line& l, double x, double y,double a)
{ l = rat_line(l.to_line().rotate(point(x, y), a), 20); }

void geowin_Rotate(rat_circle& ci, double x, double y,double a)
{ ci = rat_circle(ci.to_circle().rotate(point(x, y), a), 20); }

void geowin_Rotate(rat_polygon& p, double x, double y,double a)
{ p = rat_polygon(p.to_polygon().rotate(point(x, y), a), 20); }

void geowin_Rotate(rat_gen_polygon& p, double x, double y,double a)
{ p = rat_gen_polygon(p.to_gen_polygon().rotate(point(x, y), a), 20); }

void geowin_Rotate(d3_rat_point& p, double x, double y,double a) 
{ return; }

void geowin_Rotate(rat_rectangle& obj, double x, double y, double a)
{ return; }

#if (__LEDA__ >= 420)
void geowin_Rotate(rat_triangle& obj, double x, double y, double a)
{
  triangle t = obj.to_float().rotate(point(x,y),a);
  obj = rat_triangle(rat_point(t.point1(),20), rat_point(t.point2(),20), rat_point(t.point3(),20));
}
#endif


// *********************************************`***************************
// ******************** begin intersects **********************************

bool geowin_IntersectsBox(const rat_point& p, double x1, double y1, 
		       double x2, double y2 ,bool filled)
{ return geowin_IntersectsBox( p.to_point(), x1, y1, x2, y2, filled ); }

bool geowin_IntersectsBox(const rat_segment& s,double x1, double y1, 
		       double x2, double y2,bool filled ) 
{ return geowin_IntersectsBox( s.to_segment(),x1, y1, x2, y2, filled ); }

bool geowin_IntersectsBox(const rat_ray& r, double x1, double y1, 
		       double x2, double y2,bool filled )
{ return geowin_IntersectsBox(r.to_ray(), x1, y1, x2, y2, filled ); }

bool geowin_IntersectsBox(const rat_line& l, double x1, double y1, 
		       double x2, double y2, bool filled )
{ return geowin_IntersectsBox(l.to_line(), x1, y1, x2, y2, filled ); }

bool geowin_IntersectsBox(const rat_circle& c, double x1, double y1, 
		       double x2, double y2, bool filled ) 
{ return geowin_IntersectsBox( c.to_circle(), x1, y1, x2, y2, filled ); }

bool geowin_IntersectsBox(const rat_polygon& p, double x1, double y1, 
		       double x2, double y2, bool filled )  
{ return geowin_IntersectsBox( p.to_polygon(), x1, y1, x2, y2, filled ); }

bool geowin_IntersectsBox(const rat_gen_polygon& p, double x1, double y1, 
		       double x2, double y2, bool filled )  
{ return geowin_IntersectsBox( p.to_gen_polygon(), x1, y1, x2, y2, filled ); }

#if (__LEDA__ >= 420)
bool geowin_IntersectsBox(const rat_triangle& obj, double x1,double y1,double x2, double y2,bool f)
{ return geowin_IntersectsBox(obj.to_float(),x1,y1,x2,y2,f); }
#endif

bool geowin_IntersectsBox(const rat_rectangle& obj, double x1,double y1,double x2, double y2,bool f)
{ return geowin_IntersectsBox(obj.to_rectangle(),x1,y1,x2,y2,f); }

bool geowin_IntersectsBox(const d3_rat_point& p,double x1,double y1,double x2, double y2,bool f)
{ 
  int dm = GeoWin::get_projection_mode();
  double xw, yw;
  
  switch(dm){
   case 0:   { xw=p.xcoord().to_double(); yw=p.ycoord().to_double(); break; }
   case 1:   { xw=p.xcoord().to_double(); yw=p.zcoord().to_double(); break; }
   default: { xw=p.ycoord().to_double(); yw=p.zcoord().to_double(); break; }   
  }
  
  if (x1<=xw && x2>=xw && y1<=yw && y2>=yw) return true; else return false;
}



// ********************* end  intersects **********************************
// ************************************************************************


// ************************************************************************
// ******************** begin bounding box ********************************



void geowin_BoundingBox(const rat_point& p, double& x1, double& x2, 
		double& y1, double& y2)
{ geowin_BoundingBox( p.to_point(), x1, x2, y1, y2); }

void geowin_BoundingBox(const rat_segment& s, double& x1, double& x2,
		double& y1, double& y2)
{ geowin_BoundingBox( s.to_segment(), x1, x2, y1, y2); }

void geowin_BoundingBox(const rat_ray& r, double& x1, double& x2,
		double& y1, double& y2)
{ geowin_BoundingBox( r.to_ray(), x1, x2, y1, y2); }

void geowin_BoundingBox(const rat_line& l, double& x1, double& x2,
		double& y1, double& y2)
{ geowin_BoundingBox( l.to_line(), x1, x2, y1, y2); }

void geowin_BoundingBox(const rat_circle& c, double& x1, double& x2,
		double& y1, double& y2)
{ geowin_BoundingBox( c.to_circle(), x1, x2, y1, y2); }


void geowin_BoundingBox(const rat_polygon& p, double& x1, double& x2,
		double& y1, double& y2)
{ geowin_BoundingBox( p.to_polygon(), x1, x2, y1, y2); }

void geowin_BoundingBox(const rat_gen_polygon& p, double& x1, double& x2,
		double& y1, double& y2)
{ geowin_BoundingBox( p.to_gen_polygon(), x1, x2, y1, y2); }

#if (__LEDA__ >= 420)
void geowin_BoundingBox(const rat_triangle& obj, double& x1, double& x2,double& y1, double& y2)
{  geowin_BoundingBox(obj.to_float(),x1,x2,y1,y2); }
#endif

void geowin_BoundingBox(const rat_rectangle& obj, double& x1, double& x2,double& y1, double& y2)
{  geowin_BoundingBox(obj.to_rectangle(),x1,x2,y1,y2); }

void geowin_BoundingBox(const d3_rat_point& p, double& x1, double& x2, double& y1, double& y2)
{ 
  int dm = GeoWin::get_projection_mode();

  switch(dm){  
    case 0:   { x1=p.xcoordD(); x2=p.xcoordD(); y1=p.ycoordD(); y2=p.ycoordD(); break; }
    case 1:   { x1=p.xcoordD(); x2=p.xcoordD(); y1=p.zcoordD(); y2=p.zcoordD(); break; }
    default: { x1=p.ycoordD(); x2=p.ycoordD(); y1=p.zcoordD(); y2=p.zcoordD(); break; }   
  } 
}


// ********************* end  bounding box ********************************
// ************************************************************************


void set_rational_precision(int prec)
{
  rat_point::default_precision = prec;
}

#else

void set_rational_precision(int)
{}

#endif // NO_RAT_ALGORITHMS

GEOWIN_END_NAMESPACE
