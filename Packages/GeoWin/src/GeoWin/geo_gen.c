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
// file          : src/GeoWin/geo_gen.c
// package       : GeoWin (1.2.2)
// revision      : 1.2.2
// revision_date : 30 January 2001 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ============================================================================

#include <LEDA/geowin_init.h>
#include<LEDA/d3_rat_point.h>
#include<LEDA/rat_circle.h>
#include<LEDA/random_point.h>
#include<LEDA/random_rat_point.h>
#include<LEDA/logo.h>

GEOWIN_BEGIN_NAMESPACE

void random_points(list<point>& L, int N, double xmin, double ymin,
		   double xmax, double ymax, int grid)
{
  int x_min = (int)xmin;
  int y_min = (int)ymin;
  
  if( grid )
    {
      x_min = ((int)(xmin/grid))*grid;
      y_min = ((int)(ymin/grid))*grid;
      x_min += xmin>0.0 ? grid : 0;
      y_min += ymin>0.0 ? grid : 0;
    }
  
  int gridx = grid ? (int)(xmax-x_min)/grid : (int)(xmax-xmin);
  int gridy = grid ? (int)(ymax-y_min)/grid : (int)(ymax-ymin);
  
  int dg      = grid ? grid : 1;
  
  for(int i = 0; i < N; i++)
    { 
      double x1 = x_min + rand_int(0, gridx)*dg;
      double y1 = y_min + rand_int(0, gridy)*dg;
      
      L.append(point(x1,y1));
    }
}

void disc_points(list<point>& L, int N, circle C)
{
  if (C.is_trivial()) return;
  
  list<point> act;

  random_points_in_disc(N,100,act);
  double radius = C.radius();
  double f = 100/radius;
  double xc,yc;
  double xm=C.center().xcoord(),ym=C.center().ycoord();
  point p;
  
  forall(p,act){
    xc=xm+p.xcoord()/f;
    yc=ym+p.ycoord()/f;
    L.append(point(xc,yc));
  }
}


void lattice_points(list<point>& L, int N, double xmin, double ymin,
		   double xmax, double ymax)
{
  list<point> act;
  point p;
  double scale_x = (xmax-xmin)/220;
  double scale_y = (ymax-ymin)/220;
  
  lattice_points(N,100,act);
  forall(p,act){
    L.append(point((p.xcoord()+100)*scale_x + scale_x*10+xmin,(p.ycoord()+100)*scale_y+ scale_y*10+ymin));
  }
}

void segment_points( list<point>& L, segment s, int n)
{
  double x  = n*s.xcoord1();
  double y  = n*s.ycoord1();
  double dx = s.xcoord2() - s.xcoord1();
  double dy = s.ycoord2() - s.ycoord1();
  
  for(int i=0; i < n; i++)
  { 
    point p(x/n, y/n);
    L.append(p);
    x += dx; y += dy;
  }
}

void circle_points(list<rat_point>& L, circle C, int n)
{
  point p = C.point1();
  point q = C.point2();
  point r = C.point3();
  
  rat_point rp(p), rq(q), rr(r);
  rat_circle R(rp,rq,rr);
  
  double d = (2*LEDA_PI)/n;
  double eps = 0.001;
  double a = 0;
  
  for(int i=0; i < n; i++)
    { 
      rat_point pp = R.point_on_circle(a,eps);
      L.append(pp);
      a += d;
    }
}

void circle_segments(list<segment>& LS, circle C, int n)
{
  list<rat_point> L;
  circle_points(L, C, n);

  list_item lit = L.first();
  
  while(lit && L.succ(lit))
  {
    LS.append(segment(L[lit].to_point(), L[L.succ(lit)].to_point()));
    lit = L.succ(lit);
  }
  LS.append(segment(L.tail().to_point(), L.head().to_point()));    
}

void rectangle_points(list<point>& L, rectangle R, int n)
{
  double xmin = R.xmin();
  double xmax = R.xmax();
  double ymin = R.ymin();
  double ymax = R.ymax();
  double dx = xmax - xmin;
  double dy = ymax - ymin;
  
  random_source SR;
  double wt_x,wt_y;
  int i;
  for(i=0;i<n;i++){
    SR >> wt_x; SR >> wt_y;
    L.append(point(xmin + wt_x*dx, ymin + wt_y*dy));
  }
}

void polygon_points(list<point>& L, polygon P, int n)
{
  list<segment> LS = P.edges();
  segment siter;
  forall(siter,LS){
    list<point> Lh;
    segment_points(Lh,siter,n);
    L.conc(Lh);
  }
}

static void A(int, list<point>&, double&, double&, double, double);
static void B(int, list<point>&, double&, double&, double, double);
static void C(int, list<point>&, double&, double&, double, double);
static void D(int, list<point>&, double&, double&, double, double);

static void A(int i, list<point>& L, double& x, double& y, double dx, double dy)
{ 
  if (i > 0)
    { 
      D(i-1,L,x,y,dx,dy); x-=dx; L.append(point(x,y));
      A(i-1,L,x,y,dx,dy); y-=dy; L.append(point(x,y));
      A(i-1,L,x,y,dx,dy); x+=dx; L.append(point(x,y));
      B(i-1,L,x,y,dx,dy);
    }
}

static void B(int i, list<point>& L, double& x, double& y, double dx, double dy)
{ 
  if (i > 0)
    {
      C(i-1,L,x,y,dx,dy); y+=dy; L.append(point(x,y));
      B(i-1,L,x,y,dx,dy); x+=dx; L.append(point(x,y));
      B(i-1,L,x,y,dx,dy); y-=dy; L.append(point(x,y));
      A(i-1,L,x,y,dx,dy);
    }
}

static void C(int i, list<point>& L, double& x, double& y, double dx, double dy)
{ 
  if (i > 0)
    { 
      B(i-1,L,x,y,dx,dy); x+=dx; L.append(point(x,y));
      C(i-1,L,x,y,dx,dy); y+=dy; L.append(point(x,y));
      C(i-1,L,x,y,dx,dy); x-=dx; L.append(point(x,y));
      D(i-1,L,x,y,dx,dy);
    }
}

static void D(int i, list<point>& L, double& x, double& y, double dx, double dy)
{ 
  if (i > 0)
    { 
      A(i-1,L,x,y,dx,dy); y-=dy; L.append(point(x,y));
      D(i-1,L,x,y,dx,dy); x-=dx; L.append(point(x,y));
      D(i-1,L,x,y,dx,dy); y+=dy; L.append(point(x,y));
      C(i-1,L,x,y,dx,dy);
    }
}

void hilbert(list<point>& L, int n, double xmin, double ymin,
	     double xmax, double ymax, int grid)
{
  int d = (1 << n);

  double dx = grid ? grid : (xmax-xmin) / (d+1) ;
  double dy = grid ? grid : (ymax-ymin) / (d+1) ;
  
  double x  = grid ? ((int)(xmin/grid))*grid + d*dx : xmin + d*dx;
  double y  = grid ? ((int)(ymin/grid))*grid + d*dy : ymin + d*dy;
  
  L.append(point(x+dx,y));
  L.append(point(x,y));

  A( n, L, x, y, dx, dy);

  L.append(point(x+dx,y));
}

// ************************************************************************
// *********************** PointScene *************************************

void geowin_generate_objects(GeoWin& gw, list<point>& L)      
{
  // get calling scene ...
  // add here special handling for namespaces ??
  geo_scene c_scene = GeoWin::get_call_scene();
  //cout << "calling scene:" << c_scene << "\n";
  
  //get added generators of the scene ...
  list<string> gen_names = gw.get_generator_names(c_scene);
  
  //cout << gen_names << "\n";

  window& w = gw.get_window();
  int r    = GEOWIN_MARGIN;
  int grid = w.get_grid_mode();

  double d = r/w.scale();
  double d1 = d + 16/w.scale();
  double x1 = w.xmin()+d;
  double y1 = w.ymin()+d;
  double x2 = w.xmax()-d;
  double y2 = w.ymax()-d1;

  panel P("Points");
  
  segment seg;
  circle  cir;
  rectangle rect;
  polygon pol;
  
  int   ps = 20;
  int   pc = 20;
  int   gen_count = 0, gen_count2 = 20;
  
  P.int_item("number*100", gen_count, 0, 19);
  P.int_item("+ n", gen_count2, 0, 99);   
  
  P.button(" RANDOM ", 1);

  P.int_item("points on segment", ps, 2, 100);
  P.button(" NEAR SEGMENT ", 2);

  P.int_item("points on circle", pc, 2, 100);
  P.button(" CIRCLE ", 3);
  P.button(" LATTICE ",4);
  P.button(" IN DISC ",5);  
  P.button(" IN RECT ",6);
  P.button(" ON POLYGON ",7);
  
  int cnt = 8;
  string siter;
  forall(siter, gen_names) P.button(siter, cnt++);
  
  P.button(" DONE ", 0);

  int but;
  list<rat_point> Lr;

  but = gw.open_panel(P);
    
  switch(but)
    {
    case 0 : break;
    case 1 :
     { random_points( L, gen_count*100+gen_count2, x1, y1,  x2,  y2, grid);
       break;
     }
    case 2 :
     {
      w >> seg;
      segment_points(L, seg, ps);
      break;
     }
    case 3 :
     {
      w >> cir;
      circle_points(Lr, cir, pc);
      rat_point p;
      forall(p, Lr) L.append(p.to_point());
      break;
     }
    case 4 :
     {
      lattice_points(L, gen_count*100+gen_count2, x1, y1,  x2,  y2);
      break;
     }
    case 5 :
     {
      w >> cir;
      disc_points(L, gen_count*100+gen_count2, cir);
      break;
     }
    case 6 :
     {
      w >> rect;
      rectangle_points(L, rect, gen_count*100+gen_count2);
      break;
     }     
    case 7 :
     {
      w >> pol;
      polygon_points(L, pol, ps);
      break;
     }
    default : //cout << but << endl; break;
     { // an additional generator was choosen ...
       gw.set_generator_number(c_scene,but-8);
       break;       
     }
  }
}


// ************************************************************************
// *********************** SegmentScene ***********************************


void geowin_generate_objects(GeoWin& gw, list<segment>& S)
{ 
  geo_scene c_scene = GeoWin::get_call_scene();
  list<string> gen_names = gw.get_generator_names(c_scene);
 
  window& w = gw.get_window();
  int r     = GEOWIN_MARGIN;
  int grid  = w.get_grid_mode();

  double d2 = r/w.scale();
  double d1 = d2 + 16/w.scale();
  double x1 = w.xmin()+d2;
  double y1 = w.ymin()+d2;
  double x2 = w.xmax()-d2;
  double y2 = w.ymax()-d1;
  
  panel P("Segments");
  
  int   pc = 20;

  circle cir; 
  rectangle rect;

  int wi = 0, d = (int)(x2-x1)/10;
  int gen_count = 100;

  P.int_item("number", gen_count, 0, 300);
  P.int_item("direction", wi, 0, 360);
  P.int_item("distance", d, 0, (int)(x2-x1));
  P.int_item("segments in circle", pc, 2, 100);
  
  P.button(" RANDOM ", 1);
  P.button(" ISO_ORIENT ", 2);
  P.button(" CIRCLE ", 3);
  P.button(" IN CIRCLE ",4);  
  P.button(" IN RECT ",5);
  P.button(" DONE ", 0);
  
  int cnt = 6;
  string siter;
  forall(siter, gen_names) P.button(siter, cnt++);  
  
  list<point> L;
  list<rat_point> rL;

  list_item lit;

  int but = gw.open_panel(P);

  switch(but)
    {
    case 0 : break;
    case 1 :
     { random_points( L, 2*gen_count, x1, y1,  x2,  y2, grid);
      lit=L.first();
      while(lit)
	{
	  S.append(segment(L[lit], L[L.succ(lit)]));
	  lit = L.succ(lit); lit = L.succ(lit);
	}
      break; 
     }
    case 2 :
     { random_points( L, gen_count, x1, y1,  x2,  y2, grid);
      point p;
      forall(p, L)
	S.append(segment(p, p.translate_by_angle(2*LEDA_PI*wi/360, d)));
      break;
     }
    case 3 :
     { w >> cir;
      circle_points(rL, cir, pc);
      lit=rL.first();
      while(lit && rL.succ(lit))
	{
	  S.append(segment(rL[lit].to_point(), 
			   rL[rL.succ(lit)].to_point()));
	  lit = rL.succ(lit);
	}
      S.append(segment(rL.tail().to_point(), rL.head().to_point()));
      break;
     }

    case 4 :
     {
      w >> cir;
      disc_points(L, gen_count*2,cir);
      list_item it;
      bool fl=true;
      forall_items(it, L){
	if (fl) { S.append(segment(L[it],L[L.cyclic_succ(it)])); fl=false; }
	else fl=true;
      }
      break;
     }
    case 5 :
     {
      w >> rect;
      rectangle_points(L, rect, gen_count*2);
      list_item it;
      bool fl=true;
      forall_items(it, L){
	if (fl) { S.append(segment(L[it],L[L.cyclic_succ(it)])); fl=false; }
	else fl=true;
      }
      break;
     }     
    default : 
     { // an additional generator was choosen ...
       gw.set_generator_number(c_scene, but-6);
       break;       
     }
    }
}

void geowin_generate_objects(GeoWin& gw, list<ray>& L)
{
  list<segment> L1;
  geowin_generate_objects(gw, L1);
  segment s;
  forall(s, L1) L.append(ray(s));
}


void geowin_generate_objects(GeoWin& gw, list<line>& L)
{
  list<segment> L1;
  geowin_generate_objects(gw, L1);
  segment s;
  forall(s, L1) L.append(line(s));
}

 

// ************************************************************************
// *********************** CircleScene ************************************

void geowin_generate_objects(GeoWin& gw, list<circle>& C) 
{   
  geo_scene c_scene = GeoWin::get_call_scene();
  list<string> gen_names = gw.get_generator_names(c_scene);
  
  window& w = gw.get_window();
  int ra     = GEOWIN_MARGIN;
  int grid  = w.get_grid_mode();
 
  double d2 = ra/w.scale();
  double d1 = d2 + 16/w.scale();
  double x1 = w.xmin()+d2;
  double y1 = w.ymin()+d2;
  double x2 = w.xmax()-d2;
  double y2 = w.ymax()-d1;
  
  panel P("Circles");
  
  int gen_count = 20;

  P.int_item("number", gen_count, 0, 100);

  P.button(" RANDOM ", 1);
  P.button(" DONE ", 0);
  
  int cnt = 2;
  string siter;
  forall(siter, gen_names) P.button(siter, cnt++);


  int but = gw.open_panel(P);

  switch(but)
    {
    case 0 : break;
    case 1 :
     { list<point> L;
      random_points( L, gen_count, x1, y1,  x2,  y2, grid);
      point p;
      forall(p, L)
	{
	  int dg  = grid ? grid : 1;
	  double minr = p.xcoord()-x1;
	  if( x2 - p.xcoord() < minr ) minr = x2 - p.xcoord();
	  if( p.ycoord() - y1 < minr ) minr = p.ycoord() - y1;
	  if( y2 - p.ycoord() < minr ) minr = y2 - p.ycoord();
	  int rr = grid ? (int)(minr/grid): (int)minr;
          rr=rr+3;
	  double r = rand_int(2, rr)*dg;
	  C.append( circle(p, r) );
	}
      break;
     }
    default : 
     { // an additional generator was choosen ...
       gw.set_generator_number(c_scene, but-2);
       break;       
     }
    }
}



// ************************************************************************
// *********************** PolygonScene ***********************************

void geowin_generate_objects(GeoWin& gw, list<polygon>& Pl) 
{  
  geo_scene c_scene = GeoWin::get_call_scene();
  list<string> gen_names = gw.get_generator_names(c_scene);

  window& w = gw.get_window();
  int r     = GEOWIN_MARGIN;
  int grid  = w.get_grid_mode();

  double d2 = r/w.scale();
  double d1 = d2 + 16/w.scale();
  double x1 = w.xmin()+d2;
  double y1 = w.ymin()+d2;
  double x2 = w.xmax()-d2;
  double y2 = w.ymax()-d1;
  
  panel P("Polygons");
  int size = 3;
  int gen_count = 3;
  
  P.int_item("size", size, 3, 100);
  P.int_item("hilbert number", gen_count, 2, 6);

  P.button(" RANDOM ", 3);
  P.button(" REGULAR ", 2);
  P.button(" HILBERT ", 1);
  P.button(" DONE ", 0);
  
  int cnt = 4;
  string siter;
  forall(siter, gen_names) P.button(siter, cnt++);
  
  
  list<point>         L;
  list<point>        L1;
  list<point>        L2;
  list<rat_point>    rL;
  
  int but = gw.open_panel(P);

  switch(but)
    {
    case 0 : break;
    case 1 :
     { hilbert( L, gen_count, x1, y1,  x2,  y2, grid);
      Pl.append( polygon(L) );
      break;
     }
    case 2 :
     { circle cir;
      w >> cir;
      circle_points(rL, cir, size);
      rat_point pr;
      forall(pr, rL) L.append(pr.to_point());
      Pl.append(polygon(L));
      break;
     }
    case 3 :
     {
      random_points( L, size, x1, y1,  x2,  y2, grid);
      L.sort();
      point p1 = L.pop();
      point p2 = L.tail();
      L1.append(p1);
      while(!L.empty())
	{
	  point p = L.pop();
	  if( orientation(p1, p2, p) > 0 )
	    L1.append(p);
	  else L2.push(p);
	}
      L1.conc(L2);
      L1.reverse();
      Pl.append(polygon(L1));
      break;
     }
    default : 
     { // an additional generator was choosen ...
       gw.set_generator_number(c_scene, but-4);
       break;       
     }
    }
}

void geowin_generate_objects(GeoWin& gw, list<gen_polygon>& Pl) 
{
  geo_scene c_scene = GeoWin::get_call_scene();
  list<string> gen_names = gw.get_generator_names(c_scene);

  window& w = gw.get_window();
  int r     = GEOWIN_MARGIN;
  int grid  = w.get_grid_mode();

  double d2 = r/w.scale();
  double d1 = d2 + 16/w.scale();
  double x1 = w.xmin()+d2;
  double y1 = w.ymin()+d2;
  double x2 = w.xmax()-d2;
  double y2 = w.ymax()-d1;
  
  panel P("Polygons");
  int size = 3;
  int gen_count = 3;
  int logo_height = 2;
  
  P.int_item("size", size, 3, 100);
  P.int_item("hilbert number", gen_count, 2, 6);
  P.int_item("logo height",logo_height,1,2);

  P.button(" RANDOM ", 3);
  P.button(" REGULAR ", 2);
  P.button(" HILBERT ", 1);
  P.button(" LOGO ",4);
  P.button(" DONE ", 0);
  
  int cnt = 5;
  string siter;
  forall(siter, gen_names) P.button(siter, cnt++);
  
  list<point>        L;
  point              p;

  gen_polygon gp;
  
  int but = gw.open_panel(P);
  switch(but)
    {
    case 0 : break;
    case 1 :
     { hilbert( L, gen_count, x1, y1,  x2,  y2, grid);
      Pl.append( gen_polygon(L) );
      break;
     }
    case 2 :
    { circle  cir;
      w >> cir;
      list<rat_point> rL;
      circle_points(rL, cir, size);
      rat_point  pr;
      forall(pr, rL) L.append(pr.to_point());
      Pl.append(gen_polygon(L));
      break;
     }
    case 3 :
     { random_points( L, size, x1, y1,  x2,  y2, grid);
      L.sort();
      point p1 = L.pop(); 
      point p2 = L.tail();
      list<point> L1;
      list<point> L2;
      L1.append(p1);
      while(!L.empty())
      { point p = L.pop();
	if( orientation(p1, p2, p) > 0 )
	   L1.append(p);
	else 
           L2.push(p);
       }
      L1.conc(L2);
      Pl.append(gen_polygon(L1));
      break;
     }
    case 4 :
     { gp= leda_logo((double)logo_height); 
      Pl.append(gp);
      break;
     }
    default : 
     { // an additional generator was choosen ...
       gw.set_generator_number(c_scene, but-5);
       break;       
     }
    }  
}

#if (__LEDA__ >= 420)
void geowin_generate_objects(GeoWin& gw, list<triangle>& L)
{
}
#endif

void geowin_generate_objects(GeoWin& gw, list<rectangle>& L)
{
}

list<d3_point> conversion(list<d3_rat_point>& L,const point& wt)
{
 list<d3_point> Lf;
 d3_rat_point p;
 d3_point pfl;
 forall(p,L) {
   pfl= p.to_d3_point().translate(wt.xcoord(),wt.ycoord(),0);
   Lf.append(pfl);
 }
 return Lf;
}


void geowin_generate_objects(GeoWin& gw, list<d3_point>& L)
{
  geo_scene c_scene = GeoWin::get_call_scene();
  list<string> gen_names = gw.get_generator_names(c_scene);

  window& w = gw.get_window();
  int r     = GEOWIN_MARGIN;
  //int grid  = w.get_grid_mode();

  double d = r/w.scale();
  double x1 = w.xmin()+d;
  double y1 = w.ymin()+d;

  panel P("D3 Points");
  
  int   number = 20, number2 =0;
  
  P.int_item("number*100", number2, 0, 19);
  P.int_item("+ n", number, 0, 99);
    
  P.button(" CUBE ", 1);

  //P.int_item("max coord", maxc, 2, 100);
  P.button(" BALL ", 2);
  P.button(" SQUARE ", 3);
  P.button(" PARABOLOID ", 4);
  P.button(" LATTICE ", 5);
  P.button(" SPHERE ", 6);
  P.button(" SEGMENT ", 7);
  
  int cnt = 8;
  string siter;
  forall(siter, gen_names) P.button(siter, cnt++);   

  P.button(" DONE ", 0);

  int but;
  list<d3_rat_point> Lr;

  but = gw.open_panel(P);
  point cpoint(x1,y1);
  circle circ;

  switch (but) {
   case 0: return;
   case 1: {
     w >> circ;
     cpoint = circ.center();
     random_points_in_cube(number+number2*100,(int)(circ.radius()),Lr); 
     break;
   }
   case 2: {
     w >> circ;
     cpoint = circ.center();
     random_points_in_ball(number+number2*100,(int)(circ.radius()),Lr); 
     break;
   }
   case 3: {
     w >> circ;
     cpoint = circ.center();
     random_points_in_square(number+number2*100,(int)(circ.radius()),Lr); 
     break;
   }
   case 4: {
     w >> circ;
     cpoint = circ.center();
     random_points_on_paraboloid(number+number2*100,(int)(circ.radius()),Lr); 
     break;
   }
   case 5: { 
     w >> circ;
     cpoint = circ.center();
     lattice_points(number+number2*100,(int)(circ.radius()),Lr); 
     break; 
   }
   case 6: {
     w >> circ;
     cpoint = circ.center();
     random_points_on_sphere(number+number2*100,(int)(circ.radius()),Lr); 
     break;
   }
   case 7: {
     w >> circ;
     cpoint = circ.center();
     random_points_on_segment(number+number2*100,(int)(circ.radius()),Lr); 
     break;
   }
   default : 
     { // an additional generator was choosen ...
       gw.set_generator_number(c_scene, but-8);
       break;       
     }   
  }

  L=conversion(Lr,cpoint);
}

#if !defined(NO_RAT_ALGORITHMS)

void rat_segment_points( list<rat_point>& L, rat_segment s, int n)
{
  rational x  = n*s.xcoord1();
  rational y  = n*s.ycoord1();
  rational dx = s.xcoord2() - s.xcoord1();
  rational dy = s.ycoord2() - s.ycoord1();
  
  for(int i=0; i < n; i++)
  { 
    rat_point p(x/n, y/n);
    L.append(p);
    x += dx; y += dy;
  }
}

void rat_polygon_points(list<rat_point>& L, rat_polygon P, int n)
{
  list<rat_segment> LS = P.edges();
  rat_segment siter;
  forall(siter,LS){
    list<rat_point> Lh;
    rat_segment_points(Lh,siter,n);
    L.conc(Lh);
  }
}

void geowin_generate_objects(GeoWin& gw, list<rat_point>& L)
{
  list<point> L1;
  geowin_generate_objects(gw, L1);
  point p;
  forall(p, L1) L.append(rat_point(p));
}

void geowin_generate_objects(GeoWin& gw, list<rat_segment>& L)
{
  list<segment> L1;
  geowin_generate_objects(gw, L1);
  segment s1;
  forall(s1, L1)  L.append(rat_segment(s1));
}

void geowin_generate_objects(GeoWin& gw, list<rat_ray>& L)
{
  list<ray> L1;
  geowin_generate_objects(gw, L1);
  ray rr;
  forall(rr, L1) L.append(rat_ray(rr));
}

void geowin_generate_objects(GeoWin& gw, list<rat_line>& L)
{
  list<line> L1;
  geowin_generate_objects(gw, L1);
  line l;
  forall(l, L1) L.append(rat_line(l));
}

void geowin_generate_objects(GeoWin& gw, list<rat_circle>& L)
{
  list<circle> L1;
  geowin_generate_objects(gw, L1);
  circle c;
  forall(c, L1) L.append(rat_circle(c));
}

void geowin_generate_objects(GeoWin& gw, list<rat_polygon>& L)
{
  list<polygon> L1;
  geowin_generate_objects(gw, L1);
  polygon p;
  forall(p, L1) L.append(rat_polygon(p));
}

void geowin_generate_objects(GeoWin& gw, list<rat_gen_polygon>& Pl) 
{  
  list<gen_polygon> L1;

  geowin_generate_objects(gw, L1);
  gen_polygon p;
  forall(p, L1) Pl.append(rat_gen_polygon(p)); 
}

#if (__LEDA__ >= 420)
void geowin_generate_objects(GeoWin& gw, list<rat_triangle>& L)
{ }
#endif

void geowin_generate_objects(GeoWin& gw, list<rat_rectangle>& L) 
{ }

void geowin_generate_objects(GeoWin& gw, list<d3_rat_point>& L) 
{ 
  list<d3_point> L1;

  geowin_generate_objects(gw, L1);
  d3_point p;
  rational x,y,z;
  forall(p, L1) {
     x=rational(p.xcoord()); y=rational(p.ycoord());
     z=rational(p.zcoord());
     L.append(d3_rat_point(x,y,z)); 
  }
}

GEOWIN_END_NAMESPACE

#endif
